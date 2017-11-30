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
)

type fmt6Query struct {
	Query string
	Matches []string
}
var fmt6Column = []string{"qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"}

func BlastOffline(filename string, out string, dbname string) {
	b := blastwrapper.NcbiBlastpCommandline{}
	b.Command = `C:\Program Files\NCBI\blast-2.7.1+\bin\blastp.exe`
	b.DB = dbname
	b.Out = out
	b.OutFmt = "6"
	b.NumThreads = 4
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

func QCFmt6Query(blastCmdChan chan string, db string) {
	for b := range blastCmdChan {
		bc := blastwrapper.BlastDBCMDCommandline{}
		bc.Command = `C:\Program Files\NCBI\blast-2.7.1+\bin\blastdbcmd.exe`
		bc.DB = db
		bc.DBType = "prot"
		bc.In = b
		err := bc.Execute()
		if err != nil {
			log.Panicln(err)
		}
	}
}

func BlastFmt6Parser(filename string, outDirectory string, db string) {
	fmt6Chan := make(chan fmt6Query)
	blastCmdChan := make(chan string)
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
					fmt6Chan <- fmt6q
					break
				}
				log.Fatalln(err)
			}
			if fmt6q.Query != r[0] {
				fmt6Chan <- fmt6q
				fmt6q.Query = r[0]
				fmt6q.Matches = []string{}
			} else {
				fmt6q.Matches = append(fmt6q.Matches, r[1])
			}
		}
		close(fmt6Chan)
	} ()
	go func() {
		for fmt6 := range fmt6Chan {
			mFilename := filepath.Join(outDirectory, strings.Replace(fmt6.Query, "|", "_", -1))
			err = CreateMatchesFile(mFilename, fmt6)
			if err != nil {
				log.Panicln(err)
			} else {
				blastCmdChan <- mFilename
			}
		}
		close(blastCmdChan)
	} ()

	QCFmt6Query(blastCmdChan, db)
}

/*
func main() {
	BlastOffline(`D:\python_projects\datahelper\ancestral_wf\glycoprotein.homosapiens.fasta`, `C:\Users\localadmin\GoglandProjects\ancestral\homosapiens.fasta.blast.tsv`, `C:\Users\localadmin\GoglandProjects\ancestral\nr_customDB`)
}*/
