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
)

var root = `C:\Users\localadmin\GoglandProjects\ancestral\result`

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
}

func CreateDBHandler(w http.ResponseWriter, r *http.Request) {
	os.MkdirAll(`C:\Users\localadmin\GoglandProjects\ancestral\result\test`, os.ModePerm)
	query := workflow.LoadQuery(`D:\python_projects\datahelper\ancestral_wf\glycoprotein.homosapiens.fasta`)
	// workflow.BlastOffline(`D:\python_projects\datahelper\ancestral_wf\glycoprotein.homosapiens.fasta`, `C:\Users\localadmin\GoglandProjects\ancestral\homosapiens.fasta.blast.tsv`, `C:\Users\localadmin\GoglandProjects\ancestral\nr_customDB`)
	s := workflow.GetSpeciesList(`D:\python_projects\datahelper\ancestral_wf\species.txt`)
	filtered := make(chan workflow.BlastMap)
	go workflow.BlastFmt6Parser(`C:\Users\localadmin\GoglandProjects\ancestral\homosapiens.fasta.blast.tsv`,`C:\Users\localadmin\GoglandProjects\ancestral\result\test`, `C:\Users\localadmin\GoglandProjects\ancestral\nr_customDB`, s, query, filtered)
	sem := make(chan bool, 6)
	wg := sync.WaitGroup{}
	for f := range filtered {
		wg.Add(1)
		sem <- true
		go func() {
			CreateJSON(f)
			workflow.ProcessAlignment(f.FileName)

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
			if detail {
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
			if detail {
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
	//http.Handle("/", r)
	log.Fatal(http.ListenAndServe(":8080", handlers.CORS(originsOk, headersOk, methodsOk)(r)))
}
