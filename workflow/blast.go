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
	"regexp"
	"fmt"
	"sync"
)

var domainPatt = `(\d+)\s>*(\d+)\s(.+)\.\s\{`
var regexDomain = regexp.MustCompile(domainPatt)

type BlastMap struct {
	IdFullMap map[string]string `json:"id"`
	OrganismMap map[string]string `json:"organism"`
	AccMap map[string]string `json:"accession"`
	FileName string `json:"filename,string"`
	MatchSourceID string
	SourceSeq blastwrapper.PrimeSeq
}

type FilterResult struct {
	Seq blastwrapper.PrimeSeq
	Organism string
}

type BlastDBCMDResult struct {
	Query string
	Seq blastwrapper.PrimeSeq
	Filename string
}

type fmt6Query struct {
	Query string
	Matches []string
}

type ConcurrentBlastMap struct {
	sync.RWMutex
	Items []BlastMap
}

func (cbm *ConcurrentBlastMap) Append(bm BlastMap) {
	cbm.Lock()
	defer cbm.Unlock()
	cbm.Items = append(cbm.Items, bm)
}

var fmt6Column = []string{"qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"}
var TabQueryColumns = []string{"Entry", "Entry name", "Protein names", "Gene names", "Organism", "Length", "Topological domain", "Sequence"}

func LoadQueryTab(filename string) map[string]blastwrapper.PrimeSeq {
	f, err := os.Open(filename)
	if err != nil {
		log.Panicln(err)
	}
	qcsv := csv.NewReader(f)
	qcsv.Comma = '\t'
	header, err := qcsv.Read()
	columns := make(map[string]int)
	q := make(chan blastwrapper.PrimeSeq)
	m := make(map[string]blastwrapper.PrimeSeq)
	for i, v := range header {
		for _, v2 := range TabQueryColumns {
			if v == v2 {
				columns[v2] = i
			}
		}
	}
	go func() {
		for {
			r, err := qcsv.Read()
			if err != nil {
				if err == io.EOF {
					break
				}
			}
			s := blastwrapper.PrimeSeq{}
			s.Id = r[columns["Entry"]]
			s.Seq = strings.TrimSpace(r[columns["Sequence"]])
			s.Species = r[columns["Organism"]]
			s.Name = r[columns["Entry name"]]
			length, err := strconv.Atoi(r[columns["Length"]])
			if err != nil {
				log.Panicln(err)
			}
			s.Length =  length
			var topDomain []blastwrapper.TopDom
			for _, v := range strings.SplitN(r[columns["Topological domain"]], ";", -1) {
				v = strings.TrimSpace(v)
				if strings.HasPrefix(v, "TOPO_DOM") {
					regres := regexDomain.FindAllStringSubmatch(v, -1)
					if regres != nil {
						result := regres[0]
						if len(result) == 4 {
							start, err := strconv.Atoi(result[1])
							if err != nil {
								log.Panicln(err)
							}
							end, err := strconv.Atoi(result[2])
							if err != nil {
								log.Panicln(err)
							}
							td := blastwrapper.TopDom{Start: start, Stop: end, Type: result[3]}
							topDomain = append(topDomain, td)
						}
					}
				}
			}
			s.TopDomain = topDomain
			qc := blastwrapper.SeqQualityControl(s, true)
			if qc == true {
				q <- s
			}
		}
		close(q)
	}()

	for r := range q {
		id := strings.SplitN(r.Id, " ", 2)
		m[strings.Replace(id[0], ">", "", 1)] = r
	}
	return m
}

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
						qc := blastwrapper.SeqQualityControl(s, true)
						if qc == true {
							q <- s
						}
					}
					break
				}
				log.Panicln(err)
			}
			if strings.ContainsAny(r, ">") {
				if s.Id != "" {
					qc := blastwrapper.SeqQualityControl(s, true)
					if qc == true {
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
		m[r.Id] = len(r.Seq)
	}
	return m
}

func BlastOffline(filename string, out string, dbname string, targetNum int) {
	b := blastwrapper.NcbiBlastpCommandline{}
	b.Command = `C:\Program Files\NCBI\blast-2.7.1+\bin\blastp.exe`
	b.DB = dbname
	b.Out = out
	b.OutFmt = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
	b.NumThreads = 7
	b.EValue = 0.001
	b.Query = filename
	b.MaxTargetSeqs = targetNum
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

func BlastFmt6Parser(filename string, outDirectory string, db string, organisms []string, queryMap map[string]blastwrapper.PrimeSeq, workPool int) ConcurrentBlastMap {
	log.Printf("Processing Blast FMT6 File %v", filename)
	fmt6Chan := make(chan fmt6Query)

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
					fmt6Chan <- fmt6Query{fmt6q.Query, fmt6q.Matches[:]}
				}

				fmt6q.Query = r[0]
				fmt6q.Matches = []string{}
			}
			pi, err := strconv.ParseFloat(r[2], 32)
			if pi > 60 {
				fmt6q.Matches = append(fmt6q.Matches, r[1])
			}
			//fmt6q.Matches = append(fmt6q.Matches, r[1])
		}
		if len(fmt6q.Matches) >0 {
			fmt6Chan <- fmt6q
		}
		close(fmt6Chan)
	} ()

	//CreateBlastDBCMDFile(compiledMatches, fmt6Chan)
	var bm ConcurrentBlastMap
	sem := make(chan bool, workPool)
	wg := sync.WaitGroup{}
	for fmt6 := range fmt6Chan {
		wg.Add(1)
		sem <- true
		go func() {
			ProcessFMT6Query(fmt6, outDirectory, db, organisms, queryMap, &bm)
			defer func() {
				<- sem
			}()
			wg.Done()
		}()

	}
	wg.Wait()

	log.Println("Finished: Parsing Blast Output.")
	//close(filtered)
	log.Println(bm.Items)
	return bm
}

func ProcessFMT6Query(fmt6 fmt6Query, outDirectory string, db string, organisms []string, queryMap map[string]blastwrapper.PrimeSeq, bm *ConcurrentBlastMap) {
	log.Printf("Processing %v", fmt6.Query)
	folderName := strings.Replace(fmt6.Query, "|", "_", -1)
	compiledMatches := filepath.Join(outDirectory, folderName)
	err := os.MkdirAll(compiledMatches, os.ModePerm)
	if err != nil {
		log.Panicln(err)
	}
	fName := filepath.Join(compiledMatches, "compiled.txt")
	err = CreateMatchesFile(fName, fmt6)
	if err != nil {
		log.Panicln(err)
	}
	outName := strings.Replace(fName, "txt", "fasta", -1)
	QCFmt6Query(fName, db, outName)
	r, err := FilterMatchFile(organisms, BlastDBCMDResult{fmt6.Query, queryMap[fmt6.Query], outName})
	if err != nil {
		log.Println(err)
	} else {
		bm.Append(r)
	}
}

func FilterMatchFile(organisms []string, query BlastDBCMDResult) (BlastMap, error) {
	outFile := strings.Replace(query.Filename, ".fasta", ".filtered.fasta", -1)
	matchFile, err := os.Open(query.Filename)
	if err != nil {
		log.Panicln(err)
	}
	defer matchFile.Close()
	buff := bufio.NewReader(matchFile)
	filtered, err := os.Create(outFile)
	if err != nil {
		log.Panicln(err)
	}
	defer filtered.Close()
	b := bufio.NewWriter(filtered)

	s := blastwrapper.PrimeSeq{}

	count := 1
	result := BlastMap{}
	result.AccMap = make(map[string]string)
	result.IdFullMap = make(map[string]string)
	result.OrganismMap = make(map[string]string)
	result.FileName = outFile
	result.SourceSeq = query.Seq
	o := append([]string(nil), organisms...)
	collect := make(map[string]bool)
	check, left, found := blastwrapper.SeqFilterOrganism(query.Seq, organisms, true, collect)
	if check {
		writeOut(count, result, FilterResult{query.Seq, found}, b)
		result.MatchSourceID = strconv.Itoa(count)
		count ++
	}
	//log.Println(left)

	o = left
	//log.Println(o)
	c := make(chan FilterResult)
	go ProcessOrganisms(buff, s, o, query.Seq.Length, c, collect)

	for p := range c {
		writeOut(count, result, p, b)
		count ++
	}
	log.Println(count)
	log.Println("Finished: Filtering Fasta Sequences.")
	b.Flush()

	if (count -1) >= 20 {
		return result, nil
	} else {
		if err == nil {
			err = fmt.Errorf("started: Not enough species (%v)", outFile)
		}
	}
	return result, err

}

func writeOut(count int, result BlastMap, p FilterResult, b *bufio.Writer) {
	numb := strconv.Itoa(count)
	result.OrganismMap[numb] = p.Organism
	id := strings.SplitN(p.Seq.Id, " ", 2)
	result.IdFullMap[numb] = p.Seq.Id
	result.AccMap[numb] = id[0]
	p.Seq.Id = numb
	/*if strings.Contains(strings.ToLower(query.Seq.Species), strings.ToLower(p.Organism)) && (len(p.Seq.Seq) == len(query.Seq.Seq)) {
		log.Println(query.Filename, "matched")
		result.MatchSourceID = numb
	}*/
	b.WriteString(p.Seq.ToString())
}
func ProcessOrganisms(buff *bufio.Reader, s blastwrapper.PrimeSeq, organisms []string, queryLength int, c chan FilterResult, collect map[string]bool) {
	bound := queryLength*20/100
	for {
		r, err := buff.ReadString('\n')

		if err != nil {
			if err == io.EOF {
				if s.Id != "" {
					qc := blastwrapper.SeqQualityControl(s, true)
					if qc == true {
						check, left, found := blastwrapper.SeqFilterOrganism(s, organisms, true, collect)
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
				qc := blastwrapper.SeqQualityControl(s, true)
				if qc == true {
					check, left, found := blastwrapper.SeqFilterOrganism(s, organisms, true, collect)
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
