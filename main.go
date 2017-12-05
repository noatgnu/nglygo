package main

import (
	"github.com/noatgnu/ancestral/workflow"
	"sync"
)

func test() (tes []string) {
	return tes
}

func main() {
	query := workflow.LoadQuery(`D:\python_projects\datahelper\ancestral_wf\glycoprotein.homosapiens.fasta`)
	// workflow.BlastOffline(`D:\python_projects\datahelper\ancestral_wf\glycoprotein.homosapiens.fasta`, `C:\Users\localadmin\GoglandProjects\ancestral\homosapiens.fasta.blast.tsv`, `C:\Users\localadmin\GoglandProjects\ancestral\nr_customDB`)
	s := workflow.GetSpeciesList(`D:\python_projects\datahelper\ancestral_wf\species.txt`)
	filtered := make(chan workflow.BlastMap)
	go workflow.BlastFmt6Parser(`C:\Users\localadmin\GoglandProjects\ancestral\homosapiens.fasta.blast.tsv`,`C:\Users\localadmin\GoglandProjects\ancestral\result`, `C:\Users\localadmin\GoglandProjects\ancestral\nr_customDB`, s, query, filtered)
	sem := make(chan bool, 6)
	wg := sync.WaitGroup{}
	for f := range filtered {

		wg.Add(1)
		sem <- true
		go func() {

			workflow.ProcessAlignment(f.FileName)

			defer func() {
				<- sem
			}()
			wg.Done()
		}()
	}
	wg.Wait()
}
