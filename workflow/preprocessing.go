package workflow

import (
	"bufio"
	"github.com/noatgnu/ancestral/blastwrapper"
	"io"
	"log"
	"os"
	"strings"
)

func FilterBlastDBByOrganisms(f *os.File, organisms []string, c chan blastwrapper.PrimeSeq) {

	buff := bufio.NewReader(f)
	s := blastwrapper.PrimeSeq{}
	count := 0
	for {
		s1, err := buff.ReadString('\n')
		if err == io.EOF {
			if s.Id != "" {
				if blastwrapper.SeqFilterOrganismNoRegex(s, organisms) {
					count ++
					log.Println(count)
					c <- s
				}
			}
			break
		}
		if err != nil {
			log.Fatalln(err)
		}

		if strings.ContainsAny(s1, ">") {

			if s.Id != "" {
				if blastwrapper.SeqFilterOrganismNoRegex(s, organisms) {
					count ++
					log.Println(count)
					c <- s
				}
			}
			s = blastwrapper.PrimeSeq{}
			s.Id = strings.TrimSpace(s1)
		} else {
			s.Seq += strings.TrimSpace(s1)
		}
	}
	log.Println("Finished filtering by organisms")
	close(c)
	defer f.Close()
}

func FilterBlastDB(filename string, organisms []string, output string) {
	log.Println(filename)

	c := make(chan blastwrapper.PrimeSeq)
	f, err := os.Open(filename)
	if err != nil {
		log.Fatalln(err)
	}

	go FilterBlastDBByOrganisms(f, organisms, c)
	outfile, err := os.Create(output)
	if err != nil {
		log.Fatalln(err)
	}

	writer := bufio.NewWriter(outfile)

	for s := range c {
		writer.WriteString(s.ToString())
		log.Println(s.ToString())
	}

	log.Println("Completed blast DB Filtering")
	writer.Flush()
	defer outfile.Close()
	defer f.Close()
}

func GetSpeciesList(filename string) []string {
	f, _ := os.Open(filename)
	defer f.Close()
	buff := bufio.NewScanner(f)
	var organisms []string
	for buff.Scan() {
		organisms = append(organisms, strings.TrimSpace(buff.Text()))
	}
	return organisms
}

func PreProcessDB(filename string, output string, speciesFilename string) {
	organisms := GetSpeciesList(speciesFilename)
	log.Println("Started")
	FilterBlastDB(filename, organisms, output)
}

func CreateCustomDB(filename string, outfileName string) {
	db := blastwrapper.MakeBlastDBCommandline{}
	db.In = filename
	db.Out = outfileName
	db.DB = "prot"
	db.Execute()
}

/*func main() {
	f := `D:\python_projects\datahelper\ancestral_wf\nr`
	fo := `D:\python_projects\datahelper\ancestral_wf\nr_filtered`
	s := `D:\python_projects\datahelper\ancestral_wf\species.txt`
	PreProcessDB(f, fo, s)
	o := `C:\Users\localadmin\GoglandProjects\ancestral\nr_customDB`
	CreateCustomDB(fo, o)
}*/