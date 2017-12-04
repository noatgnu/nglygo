package workflow

import (
	"github.com/noatgnu/ancestral/clustalowrapper"
	"log"
	"github.com/noatgnu/ancestral/phymlwrapper"
	"strings"
)

func CreateAlignment(in string, out string) {
	c := clustalowrapper.ClustalOCommandline{}
	c.Command = `D:\clustal-omega-1.2.2-win64\clustalo.exe`
	c.In = in
	c.OutFmt = "--outfmt=phylip"
	c.Out = out
	err := c.Execute()
	if err != nil {
		log.Panicln(err)
	}
}

func CreateTreeML(in string) {
	c := phymlwrapper.PhyMLCommandline{}
	c.Command = `D:\PhyML-3.1\PhyML-3.1_win32.exe`
	c.In = in
	err := c.Execute()
	if err != nil {
		log.Panicln(err)
	}
}

func ProcessAlignment(filename string) {
	log.Printf("Started: Phylogeny Construction (%v)", filename)
	alignmentFile := strings.Replace(filename, ".filtered.fasta", ".phy", -1)
	CreateAlignment(filename, alignmentFile)
	CreateTreeML(alignmentFile)

}