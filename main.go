package main

import (
	"github.com/noatgnu/ancestral/workflow"
)

func test() (tes []string) {
	return tes
}

func main() {
	workflow.BlastOffline(`D:\python_projects\datahelper\ancestral_wf\glycoprotein.homosapiens.fasta`, `C:\Users\localadmin\GoglandProjects\ancestral\homosapiens.fasta.blast.tsv`, `C:\Users\localadmin\GoglandProjects\ancestral\nr_customDB`)
}
