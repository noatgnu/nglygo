package workflow

import (
	"github.com/noatgnu/ancestral/blastwrapper"
	"log"
	"os"
	"encoding/csv"
	"io"
	"bufio"
	"strings"
	"path/filepath"
	"strconv"
)

type BlastMap struct {
	IdFullMap map[string]string `json:"id"`
	OrganismMap map[string]string `json:"organism"`
	AccMap map[string]string `json:"accession"`
	FileName string `json:"filename,string"`
}

type FilterResult struct {
	Seq blastwrapper.PrimeSeq
	Organism string
}

type BlastDBCMDResult struct {
	Query string
	QueryLength int
	Filename string
}

type fmt6Query struct {
	Query string
	Matches []string
}
var fmt6Column = []string{"qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"}

func LoadQuery(filename string) map[string]int {
	f, err := os.Open(filename)
	if err != nil {
		log.Panicln(err)
	}
	buff := bufio.NewReader(f)
	q := make(chan blastwrapper.PrimeSeq)
	s := blastwrapper.PrimeSeq{}
	m := make(map[string]int)
	go func() {
		for {
			r, err := buff.ReadString('\n')
			if err != nil {
				if err == io.EOF {
					if s.Id != "" {
						if blastwrapper.SeqQualityControl(s, true) {
							q <- s
						}
					}
					break
				}
				log.Panicln(err)
			}
			if strings.ContainsAny(r, ">") {
				if s.Id != "" {
					if blastwrapper.SeqQualityControl(s, true) {
						q <- s
					}
				}
				s = blastwrapper.PrimeSeq{}
				s.Id = strings.TrimSpace(r)
			} else {
				s.Seq += strings.TrimSpace(r)
			}
		}
		close(q)
	}()
	for r := range q {
		id := strings.SplitN(r.Id, " ", 2)
		m[strings.Replace(id[0], ">", "", 1)] = len(r.Seq)
	}
	return m
}

func BlastOffline(filename string, out string, dbname string) {
	b := blastwrapper.NcbiBlastpCommandline{}
	b.Command = `C:\Program Files\NCBI\blast-2.7.1+\bin\blastp.exe`
	b.DB = dbname
	b.Out = out
	b.OutFmt = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
	b.NumThreads = 7
	b.EValue = 0.001
	b.Query = filename
	err := b.Execute()
	if err != nil {
		log.Panicln(err)
	}
	log.Println("Finished.")
}

func CreateMatchesFile(filename string, fmt6q fmt6Query) error {
	f, err := os.Create(filename)
	if err != nil {
		log.Panicln(err)
	}
	defer f.Close()
	writer := bufio.NewWriter(f)
	for _, v := range fmt6q.Matches {
		writer.WriteString(v+"\n")
	}
	defer writer.Flush()
	return err
}

func CreateBlastDBCMDFile(filename string, c chan fmt6Query) {
	outfile, err := os.Create(filename)
	if err != nil {
		log.Panicln(err)
	}
	defer outfile.Close()
	writer := bufio.NewWriter(outfile)
	defer writer.Flush()
	for fmt6 := range c {
		writer.WriteString(strings.Join(fmt6.Matches, "\n")+"\n")
	}

}

func QCFmt6Query(filename string, db string, outFilename string) {
	bc := blastwrapper.BlastDBCMDCommandline{}
	bc.Command = `C:\Program Files\NCBI\blast-2.7.1+\bin\blastdbcmd.exe`
	bc.DB = db
	bc.DBType = "prot"
	bc.In = filename
	bc.Out = outFilename
	err := bc.Execute()
	if err != nil {
		log.Panicln(err)
	}
}

func BlastFmt6Parser(filename string, outDirectory string, db string, organisms []string, queryMap map[string]int, filtered chan BlastMap) {
	fmt6Chan := make(chan fmt6Query)
	BlastDBCMDChan := make(chan BlastDBCMDResult)
	f, err := os.Open(filename)
	if err != nil {
		log.Panicln(err)
	}
	defer f.Close()

	c := csv.NewReader(f)
	c.Comma = '\t'
	go func() {
		var fmt6q fmt6Query
		for {
			r, err := c.Read()
			if err != nil {
				if err == io.EOF {
					if fmt6q.Query != "" {
						fmt6Chan <- fmt6q
					}
					break
				}
				log.Fatalln(err)
			}

			if fmt6q.Query != r[0] {
				if fmt6q.Query != "" {
					fmt6Chan <- fmt6q
				}
				fmt6q.Query = r[0]
				fmt6q.Matches = []string{}
			}
			fmt6q.Matches = append(fmt6q.Matches, r[1])
		}
		close(fmt6Chan)
	} ()

	//CreateBlastDBCMDFile(compiledMatches, fmt6Chan)

	go func () {
		for fmt6 := range fmt6Chan {
			folderName := strings.Replace(fmt6.Query, "|", "_", -1)
			compiledMatches := filepath.Join(outDirectory, folderName)
			os.MkdirAll(compiledMatches, os.ModePerm)
			fName := filepath.Join(compiledMatches, "compiled.txt")
			err = CreateMatchesFile(fName, fmt6)
			if err != nil {
				log.Panicln(err)
			} else {
				outName := strings.Replace(fName, "txt", "fasta", -1)
				QCFmt6Query(fName, db, outName)
				BlastDBCMDChan <- BlastDBCMDResult{fmt6.Query, queryMap[fmt6.Query], outName}
			}
		}
	} ()

	for bcmd := range BlastDBCMDChan {
		FilterMatchFile(organisms, bcmd, filtered)
	}
	close(filtered)

}

func FilterMatchFile(organisms []string, query BlastDBCMDResult, filteredChan chan BlastMap) {
	o := append([]string(nil), organisms...)
	outFile := strings.Replace(query.Filename, ".fasta", ".filtered.fasta", -1)
	matchFile, err := os.Open(query.Filename)
	if err != nil {
		log.Panicln(err)
	}
	defer matchFile.Close()
	buff := bufio.NewReader(matchFile)
	s := blastwrapper.PrimeSeq{}
	c := make(chan FilterResult)
	go ProcessOrganisms(buff, s, o, query.QueryLength, c)

	filtered, err := os.Create(outFile)
	if err != nil {
		log.Panicln(err)
	}
	defer filtered.Close()
	b := bufio.NewWriter(filtered)
	count := 1
	result := BlastMap{}
	result.AccMap = make(map[string]string)
	result.IdFullMap = make(map[string]string)
	result.OrganismMap = make(map[string]string)
	result.FileName = outFile
	for p := range c {
		numb := strconv.Itoa(count)
		result.OrganismMap[numb] = p.Organism
		id := strings.SplitN(p.Seq.Id, " ", 2)
		result.IdFullMap[numb] = p.Seq.Id
		result.AccMap[numb] = id[0]
		p.Seq.Id = ">"+numb
		b.WriteString(p.Seq.ToString())
		count++
	}
	if count >= 20 {
		filteredChan <- result
	} else {
		log.Printf("Started: Not enough species (%v)", outFile)
	}
	defer b.Flush()

}
func ProcessOrganisms(buff *bufio.Reader, s blastwrapper.PrimeSeq, organisms []string, queryLength int, c chan FilterResult) {
	bound := queryLength*20/100
	for {
		r, err := buff.ReadString('\n')
		if err != nil {
			if err == io.EOF {
				if s.Id != "" {
					if blastwrapper.SeqQualityControl(s, true) {
						check, left, found := blastwrapper.SeqFilterOrganism(s, organisms, true)
						if check {
							l := len(s.Seq)
							if (queryLength - bound) <= l && l <= (queryLength + bound) {
								c <- FilterResult{s, found}
								organisms = left
							}
						}
					}
				}
				break
			}
			log.Panicln(err)
		}
		if strings.ContainsAny(r, ">") {
			if s.Id != "" {
				if blastwrapper.SeqQualityControl(s, true) {
					check, left, found := blastwrapper.SeqFilterOrganism(s, organisms, true)
					if check {
						l := len(s.Seq)
						if (queryLength - bound) < l && l < (queryLength + bound) {
							c <- FilterResult{s, found}
							organisms = left
						}
					}
				}
			}
			s = blastwrapper.PrimeSeq{}
			s.Id = strings.TrimSpace(r)
		} else {
			s.Seq += strings.TrimSpace(r)
		}
	}
	close(c)
}

/*
func main() {
	BlastOffline(`D:\python_projects\datahelper\ancestral_wf\glycoprotein.homosapiens.fasta`, `C:\Users\localadmin\GoglandProjects\ancestral\homosapiens.fasta.blast.tsv`, `C:\Users\localadmin\GoglandProjects\ancestral\nr_customDB`)
}*/
