package main

import (
	"github.com/noatgnu/ancestral/workflow"
	"github.com/gorilla/mux"
	"sync"
	"os"
	"strings"
	"log"
	"encoding/json"
	"net/http"
	"io/ioutil"
	"path/filepath"
	"bufio"
	"io"
	"github.com/gorilla/handlers"
	"encoding/csv"
	"github.com/noatgnu/ancestral/blastwrapper"
	"github.com/noatgnu/ancestral/alignio"
	"strconv"
	"sort"
)

var root = `D:\GoProject\ancestral\result`

type Response struct {
	DB []CoreDB   `json:"db"`
	Query []Query `json:"query"`
}

type Query struct {
	ID string `json:"id"`
	Tree string `json:"tree"`
	HasTree bool `json:"hasTree"`
	TreeFile string `json:"treeFile"`
	HasAlign bool `json:"hasAlign"`
	AlignFile string `json:"alignFile"`
	Description string `json:"description"`
	BlastMap workflow.BlastMap `json:"blastMap"`
	HasMotifAnalysis bool `json:"hasMotif"`
	MotifAnalysisFile string `json:"motifFile"`
	MotifAnalysis []workflow.BranchMotif `json:"motif"`
	HasTopDom bool `json:"hasTopDom"`
	TopDomFile string `json:"topDomFile"`
	TopDomDF []blastwrapper.TopDom `json:"topDomDF"`
	Alignment alignio.Alignment `json:"alignment"`
	AlignmentDF []workflow.AlignmentData `json:"alignmentDF"`
}

type CoreDB struct {
	DBName string `json:"dbName"`
	Results []string `json:"results"`
}

func CreateJSON(blastMap workflow.BlastMap) {
	f, err := os.Create(strings.Replace(blastMap.FileName, ".fasta", ".json", -1))
	if err != nil {
		log.Panicln(err)
	}
	jw := json.NewEncoder(f)
	jw.Encode(blastMap)
	f.Close()
}

func CreateDBHandler(w http.ResponseWriter, r *http.Request) {
	os.MkdirAll(`D:\GoProject\ancestral\result\test2`, os.ModePerm)
	query := workflow.LoadQueryTab(`D:\python_projects\datahelper\ancestral_wf\uniprot-glycoprotein.tab`)
	// workflow.BlastOffline(`D:\python_projects\datahelper\ancestral_wf\glycoprotein.homosapiens.fasta`, `C:\Users\localadmin\GoglandProjects\ancestral\homosapiens.fasta.blast.tsv`, `C:\Users\localadmin\GoglandProjects\ancestral\nr_customDB`)
	s := workflow.GetSpeciesList(`D:\python_projects\datahelper\ancestral_wf\species.txt`)
	filtered := make(chan workflow.BlastMap)
	go workflow.BlastFmt6Parser(`D:\GoProject\ancestral\homosapiens.fasta.blast.tsv`,`D:\GoProject\ancestral\result\test2`, `D:\GoProject\ancestral\nr_customDB`, s, query, filtered)
	sem := make(chan bool, 6)
	wg := sync.WaitGroup{}
	for f := range filtered {
		CreateJSON(f)
		if f.MatchSourceID != "" {
			log.Println(f.FileName, "has match")
		} else {
			log.Println(f.FileName, "no match")
		}
		wg.Add(1)
		sem <- true
		go func() {
			workflow.ProcessAlignment(f)
			defer func() {
				<- sem
			}()
			wg.Done()
		}()
	}
	wg.Wait()
}

func HomeHandler(w http.ResponseWriter, r *http.Request) {

}

func MakeBlastDBHandler(w http.ResponseWriter, r *http.Request) {
	if r.Method == "GET" {
		f := `D:\python_projects\datahelper\ancestral_wf\nr`
		fo := `D:\python_projects\datahelper\ancestral_wf\nr_filtered`
		s := `D:\python_projects\datahelper\ancestral_wf\species.txt`
		workflow.PreProcessDB(f, fo, s)
		o := `D:\GoProject\ancestral\nr_customDB`
		workflow.CreateCustomDB(fo, o)
	}
}

func BlastPHandler(w http.ResponseWriter, r *http.Request){
	f, err := os.Open(`D:\python_projects\datahelper\ancestral_wf\uniprot-glycoprotein.tab`)
	if err != nil {
		log.Panicln(err)
	}
	qcsv := csv.NewReader(f)
	qcsv.Comma = '\t'
	header, err := qcsv.Read()
	if err != nil {
		log.Panicln(err)
	}
	columns := make(map[string]int)
	q := make(chan blastwrapper.PrimeSeq)
	for i, v := range header {
		for _, v2 := range workflow.TabQueryColumns {
			if v == v2 {
				columns[v2] = i
			}
		}
	}
	go func(){
		for {
			r, err := qcsv.Read()

			if err != nil {
				if err == io.EOF {
					break
				}
				log.Fatalln(err)
			}
			s := blastwrapper.PrimeSeq{}
			s.Id = r[columns["Entry"]]
			s.Seq = r[columns["Sequence"]]
			qc := blastwrapper.SeqQualityControl(s, true)
			if qc == true {
				q <- s
			}
		}
		defer close(q)
	}()

	o, err := os.Create(`D:\GoProject\ancestral\uniprot-homosapiens.fasta`)
	if err != nil {
		log.Panicln(err)
	}
	writer := bufio.NewWriter(o)
	count := 0
	log.Println("Parsing tabulated file.")
	for s := range q {
		count++
		log.Println(count)
		writer.WriteString(s.ToString())
	}
	log.Println("Finished parsing tabulated file.")
	writer.Flush()
	o.Close()

	if r.Method == "GET" {
		workflow.BlastOffline(`D:\GoProject\ancestral\uniprot-homosapiens.fasta`, `D:\GoProject\ancestral\homosapiens.fasta.blast.tsv`, `D:\GoProject\ancestral\nr_customDB`)
	}
}

func GetDBListHandler(w http.ResponseWriter, r *http.Request) {
	if r.Method == "GET" {
		var response Response
		f, err := ioutil.ReadDir(root)
		if err != nil {
			log.Panicln(err)
		}
		for _, v := range f {
			response.DB = append(response.DB, CoreDB{v.Name(), []string{}})
		}

		j := json.NewEncoder(w)
		j.Encode(response)
	}
}

func GetQueryHandler(w http.ResponseWriter, r *http.Request) {
	vars := mux.Vars(r)
	db := vars["db"]

	id := strings.Split(vars["id"],";")

	if r.Method == "GET" {
		var response Response
		for _, v := range id {
			if v != "" {
				query := GetQuery(v, db, true)
				response.Query = append(response.Query, query)
			}
		}

		j := json.NewEncoder(w)
		j.Encode(response)
	}
}

func GetQuery (id string, db string, detail bool) Query {
	var query Query
	query.ID = id
	f, err := filepath.Glob(filepath.Join(root, db, id, "*"))
	if err != nil {
		log.Panicln(err)
	}
	for _, v := range f {
		if strings.HasSuffix(v, `compiled.filtered.json`) {
			if detail == true {
				fj, err := os.Open(v)
				if err != nil {
					log.Panicln(err)
				}
				var m workflow.BlastMap
				j := json.NewDecoder(fj)
				j.Decode(&m)
				query.BlastMap = m
				fj.Close()
			}
		} else if strings.HasSuffix(v, `.reconstructed.tree.txt`) {
			query.TreeFile = v
			query.HasTree = true
			if detail == true{
				fj, err := os.Open(v)
				if err != nil {
					log.Panicln(err)
				}
				rd := bufio.NewReader(fj)
				for {
					r, err := rd.ReadString('\n')

					if err == io.EOF {
						break
					}
					query.Tree += r
				}
				//log.Println(query.Tree)
				fj.Close()
			}
		} else if strings.HasSuffix(v, `compiled.reconstructed.phy`) {
			query.AlignFile = v
			query.HasAlign = true
			motifFile := strings.Replace(v,".reconstructed.phy", ".reconstructed.motifs.txt", -1)
			hasMotif := false
			if _, err := os.Stat(motifFile); err == nil {
				hasMotif = true
			}
			topDomPath := strings.Replace(v,".reconstructed.phy", ".reconstructed.topdom.txt", -1)
			if _, err := os.Stat(topDomPath); err == nil {
				query.HasTopDom = true
			}
			if detail == true{
				query.Alignment = alignio.ReadPhylip(v)
				if hasMotif == true {
					o, err := workflow.ReadTabulatedFile(motifFile)
					if err != nil {
						log.Panicln(err)
					}
					if o != nil {
						query.Alignment.ConserveMap = workflow.ConserveCount(o)
					}
				}
				if query.HasTopDom == true {
					top, err := workflow.ReadTabulatedFile(topDomPath)
					if err != nil {
						log.Panicln(err)
					}
					if top != nil {
						query.HasTopDom = true
						query.TopDomFile = topDomPath
						query.TopDomDF = workflow.TopDom(top)
					}
				}

				totalrow := query.Alignment.Length*len(query.Alignment.Alignment)
				result := make([]workflow.AlignmentData, totalrow)
				rowNumb := 0
				yCoord := 0
				min := 0
				max := 0
				val := 0
				var seqidarray alignio.BySeqId

				for k := range query.Alignment.Alignment {
					seqidarray = append(seqidarray, alignio.SeqD3Coordinate{SeqId: k})
				}
				sort.Sort(seqidarray)

				for i, v2 := range seqidarray {
					seqid, err := strconv.Atoi(v2.SeqId)
					if err != nil {
						log.Panicln(err)
					}
					for pos := 0; pos < query.Alignment.Length; pos ++ {
						extra := make(map[string]string)
						for _, td := range query.TopDomDF {
							if td.Start <= (pos+1) && td.Stop >= (pos+1) {
								extra["topDomType"] = td.Type
								break
							}
						}
						result[rowNumb] = workflow.AlignmentData{SeqId: seqid, AA:  query.Alignment.Alignment[v2.SeqId][pos:pos+1], Pos: pos+1, YCoord: yCoord, Value: val, Extra:extra}
						rowNumb ++
						if val > max {
							max = val
						}
						if val < min {
							min = val
						}

					}
					seqidarray[i].YCoord = yCoord
					yCoord += 1
				}

				query.AlignmentDF = result[:]
				query.Alignment.SeqIdArray = seqidarray
			}
		} else if strings.HasSuffix(v, `compiled.motif.analysis.txt`) {
			query.MotifAnalysisFile = v
			query.HasMotifAnalysis = true
			if detail == true{
				query.MotifAnalysis = workflow.ReadMotifAnalysis(v)
			}
		}
	}
	return query
}

func GetAllQueryHandler(w http.ResponseWriter, r *http.Request) {
	vars := mux.Vars(r)
	db := vars["db"]
	var response Response
	if r.Method == "GET" {
		f, err := ioutil.ReadDir(filepath.Join(root, db))
		if err != nil {
			log.Panicln(err)
		}
		cdb := CoreDB{db, []string{}}
		for _, v := range f {
			var query Query
			query = GetQuery(v.Name(), db, false)
			response.Query = append(response.Query, query)
			cdb.Results = append(cdb.Results, v.Name())
		}
		response.DB = append(response.DB, cdb)
		j := json.NewEncoder(w)
		j.Encode(response)
	}
}

func init() {

}

func main() {
	headersOk := handlers.AllowedHeaders([]string{"X-Requested-With"})
	originsOk := handlers.AllowedOrigins([]string{"*"})
	methodsOk := handlers.AllowedMethods([]string{"GET", "HEAD", "POST", "PUT", "OPTIONS"})
	r := mux.NewRouter()
	r.HandleFunc("/", HomeHandler)
	r.HandleFunc(`/db/`, GetDBListHandler)
	r.HandleFunc(`/db/{db}`, GetAllQueryHandler)
	r.HandleFunc(`/query/{db}/{id}`, GetQueryHandler)
	r.HandleFunc(`/create`, CreateDBHandler)
	r.HandleFunc(`/makeblastdb`, MakeBlastDBHandler)
	r.HandleFunc(`/blastp`, BlastPHandler)
	//http.Handle("/", r)
	log.Fatal(http.ListenAndServe(":8080", handlers.CORS(originsOk, headersOk, methodsOk)(r)))
}
