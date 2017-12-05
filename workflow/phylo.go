package workflow

import (
	"github.com/noatgnu/ancestral/clustalowrapper"
	"log"
	"github.com/noatgnu/ancestral/phymlwrapper"
	"strings"
	"github.com/noatgnu/ancestral/codemlwrapper"
)

func CreateAlignment(in string, out string) {
	c := clustalowrapper.ClustalOCommandline{}
	c.Command = `D:\clustal-omega-1.2.2-win64\clustalo.exe`
	c.In = in
	c.OutFmt = "--outfmt=phylip"
	c.Out = out+"_inter"
	err := c.Execute()
	if err != nil {
		log.Panicln(err)
	}
	c.ConvertSequential(out+"_inter", out)
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

func ASR(seq string, tree string, out string) {
	a := codemlwrapper.CodeMLCommandline{}
	a.Command = `D:\paml4.9e\bin\codeml.exe`
	a.SeqFile = seq
	a.TreeFile = tree
	a.OutFile = out
	a.Noisy = 9
	a.Verbose = 2
	a.Runmode = 0
	a.SeqType = 2
	a.Clock = 0
	a.AADist = 0
	a.AARateFile = `D:\paml4.9e\dat\wag.dat`
	a.Model = 2
	a.ICode = 0
	a.Mgene = 0
	a.FixAlpha = 0
	a.Alpha = "0.5"
	a.MAlpha = 1
	a.NCatG = 4
	a.GetSE = 0
	a.RateAncestor = 1
	a.SmallDiff = "5e-07"
	a.CleanData = 0
	a.Method = 1
	a.BuildCtl(strings.Replace(seq, ".phy", ".ctl", -1))
	err := a.Execute()
	if err != nil {
		log.Panicln(err)
	}
}

func ProcessAlignment(filename string) {
	log.Printf("Started: Phylogeny Construction (%v)", filename)
	alignmentFile := strings.Replace(filename, ".filtered.fasta", ".phy", -1)
	CreateAlignment(filename, alignmentFile)
	CreateTreeML(alignmentFile)
	ASR(alignmentFile,alignmentFile+"_phyml_tree.txt", strings.Replace(alignmentFile, ".phy", ".asr", -1))
}