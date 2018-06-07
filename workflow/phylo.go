package workflow

import (
	"github.com/noatgnu/ancestral/clustalowrapper"
	"log"
	"github.com/noatgnu/ancestral/phymlwrapper"
	"strings"
	"github.com/noatgnu/ancestral/codemlwrapper"
	"github.com/noatgnu/ancestral/pythoncmdwrapper"
	"github.com/noatgnu/ancestral/alignio"
	"github.com/noatgnu/ancestral/nglycan"
	"os"
	"bufio"
	"fmt"
	"strconv"
	"regexp"
	"encoding/csv"
	"io"
	"github.com/noatgnu/ancestral/blastwrapper"
	"github.com/noatgnu/ancestral/fasttreecmdwrapper"
	"github.com/noatgnu/ancestral/raxmlwrapper"
	"path/filepath"
)

const nRegexPat = `N`
var nRegex = regexp.MustCompile(nRegexPat)
const gappedRegexPat = `-*\w+`
var gappedRegex = regexp.MustCompile(gappedRegexPat)
const bootstrapedTree = `\)[\d|\.]+\:`
var treeRegex = regexp.MustCompile(bootstrapedTree)


type BranchMotif struct {
	Origin string `json:"origin"`
	Target string `json:"target"`
	Source string `json:"source"`
	OriginStart string `json:"originStart"`
	OriginEnd string `json:"originEnd"`
	OriginSeq string `json:"originSeq"`
	TargetStart string `json:"targetStart"`
	TargetEnd string `json:"targetEnd"`
	TargetSeq string `json:"targetSeq"`
	ExpandedOriginStart string `json:"expOriginStart"`
	ExpandedOriginEnd string `json:"expOriginEnd"`
	ExpandedOriginSeq string `json:"expOriginSeq"`
	ExpandedTargetStart string `json:"expTargetStart"`
	ExpandedTargetEnd string `json:"expTargetEnd"`
	ExpandedTargetSeq string `json:"expTargetSeq"`
	Status string `json:"status"`
}

type AlignmentData struct {
	SeqId int `json:"seqid"`
	AA string `json:"aa"`
	Pos int	`json:"pos"`
	Value int `json:"value"`
	YCoord int `json:"yCoord"`
	Extra map[string]string `json:"extra"`
}

type MotifConservity struct {
	Pos int `json:"pos"`
	Conserve int `json:"conserve"`
	Count int `json:"count"`
}

func TopDom(inchan chan map[string]string) []blastwrapper.TopDom {
	var topDom []blastwrapper.TopDom
	for m := range inchan {

		start, err := strconv.Atoi(m["Remapped_Start"])
		if err != nil {
			log.Println(m)
			log.Panicln(err)
		}

		end, err := strconv.Atoi(m["Remapped_End"])
		if err != nil {
			log.Panicln(err)
		}
		topDom = append(topDom, blastwrapper.TopDom{Start: start, Stop: end, Type: m["Type"]})
	}
	return topDom
}

func ReadTabulatedFile(filename string) (chan map[string]string, error) {
	infile, err := os.Open(filename)
	if err != nil {
		log.Panicln(err)
	}
	o := make(chan map[string]string)
	csvreader := csv.NewReader(infile)
	csvreader.Comma = '\t'
	header, err := csvreader.Read()
	if err != nil {
		return nil, err
	}

	go func() {
		for {
			row, err := csvreader.Read()
			if err != nil {
				if err == io.EOF {
					break
				}
				log.Panicln(err)
			}
			m := make(map[string]string)
			for i, v := range row {
				m[header[i]] = v
			}
			o <- m
		}
		close(o)
		infile.Close()
	}()
	return o, nil
}

func ConserveCount(inchan chan map[string]string) map[int]int {
	conserveMap := make(map[int]int)
	for m := range inchan {
		pos, err := strconv.Atoi(m["Start_Position_With_Gap"])
		if err != nil {
			log.Panicln(err)
		}
		if _, ok := conserveMap[pos]; !ok {
			conserveMap[pos] = 0
		}
		conserveMap[pos] ++
	}
	return conserveMap
}


func CreateAlignment(in string, out string) {
	c := clustalowrapper.ClustalOCommandline{}
	c.Command = `D:\clustal-omega-1.2.2-win64\clustalo.exe`
	c.In = in

	if strings.HasSuffix(out, "phy"){
		c.OutFmt = "--outfmt=phylip"
		c.Out = out+"_inter"
	} else {
		c.OutFmt = "--outfmt=fasta"
		c.Out = out
	}
	err := c.Execute()
	if err != nil {
		log.Panicln(err)
	}

	c.ConvertAlignment(out+"_inter", out, "phylip")
	c.ConvertAlignment(out, out+".fasta", "fasta")

}

func ReadAlignmentPhylip(filename string, b BlastMap) (map[string][][]int, map[string][][]int, alignio.Alignment, map[string]map[int][]int) {
	a := alignio.ReadPhylip(filename)
	f, err := os.Create(strings.Replace(filename, ".phy", ".motifs.txt", -1))
	if err != nil {
		log.Panicln(err)
	}
	defer f.Close()
	writer := bufio.NewWriter(f)
	writer.WriteString("Entry\tStart_Position\tEnd_Position\tMotif\tStart_Position_With_Gap\tEnd_Position_With_Gap\tMotif_With_Gap\tTopological_Domain\n")

	motifMap := make(map[string][][]int)
	alignedMotifMap := make(map[string][][]int)
	nResidueAlignedMotifMap := make(map[string]map[int][]int)
	var reMap []blastwrapper.TopDom
	if b.SourceSeq.TopDomain != nil {
		if b.MatchSourceID != "" {
			reMap = RemapTopDom(a.Alignment[b.MatchSourceID], b.SourceSeq.TopDomain, strings.Replace(filename, ".phy", ".topdom.txt", -1))
		}
	}
	for k, v := range a.Alignment {
		nResidueAlignedMotifMap[k] = make(map[int][]int)
		gap := nglycan.MotifParseStringWithGap(v, -1)
		noGapSeq := strings.Replace(v, "-", "", -1)
		noGap := nglycan.MotifParseStringWithGap(noGapSeq, -1)
		motifMap[k] = noGap

		alignedMotifMap[k] = gap
		for i := range gap {
			nResidueAlignedMotifMap[k][gap[i][0]] = gap[i]
			topdomType := ""
			if reMap != nil {
				for _, rm := range reMap {
					if (rm.Start <= gap[i][0]) && (rm.Stop > gap[i][0]) {
						topdomType = rm.Type
						break
					}
				}
			}
			writer.WriteString(fmt.Sprintf("%v\t%v\t%v\t%v\t%v\t%v\t%v\t%v\n", k, noGap[i][0]+1, noGap[i][1], noGapSeq[noGap[i][0]:noGap[i][1]], gap[i][0]+1, gap[i][1], v[gap[i][0]:gap[i][1]], topdomType))
		}
	}
	writer.Flush()
	return motifMap, alignedMotifMap, a, nResidueAlignedMotifMap
}

func CreateTreeML(in string, bootstrap int, speed string) {
	switch speed {
	case "fast":
		c := fasttreecmdwrapper.FastMLCommandline{}
		c.Command = `D:\GoProject\ancestral\FastTree.exe`
		c.In = in+".fasta"
		c.Out = in+"_phyml_tree.txt"
		log.Println(c.Out)
		err := c.Execute()
		if err != nil {
			log.Panicln(err)
		}
	case "moderate":
		c := raxmlwrapper.RaxMLCommandline{}
		c.Command = `D:\GoProject\ancestral\FastTree.exe`
		c.In = in
		c.Out = in+"_phyml_tree.txt"
		log.Println(c.Out)
		err := c.Execute()
		if err != nil {
			log.Panicln(err)
		}
	case "slow":
		c := phymlwrapper.PhyMLCommandline{}
		c.Command = `D:\PhyML-3.1\PhyML-3.1_win32.exe`
		c.In = in
		c.Bootstrap = bootstrap
		err := c.Execute()
		if err != nil {
			log.Panicln(err)
		}
	}
}

func ASR(seq string, tree string, out string, bm BlastMap) {
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
	a.Verbose = 2
	a.NData = 1
	a.BuildCtl(strings.Replace(seq, ".phy", ".ctl", -1))
	err := a.Execute()
	if err != nil {
		log.Panicln(err)
	}
	dir, _ := filepath.Split(a.SeqFile)
	branches := a.ReadSupplemental(filepath.Join(dir, "rst"))
	MotifAnalysis(seq, branches, 0, bm)
	CombineTree(a)
}

func ReadMotifAnalysis(filename string) (ArrayMotifA []BranchMotif) {
	f, err := os.Open(filename)
	if err != nil {
		log.Panicln(err)
	}
	reader := csv.NewReader(f)
	reader.Comma = '\t'
	var header []string
	for {
		r, err := reader.Read()
		if err != nil {
			if err == io.EOF {
				break
			}
		}
		if header != nil {
			ArrayMotifA = append(ArrayMotifA, BranchMotif{r[0], r[1], r[2],r[3], r[4], r[5],r[6],r[7],r[8],r[9], r[10],r[11],r[12],r[13],r[14], r[15]})
		} else {
			header = r
		}
	}
	return ArrayMotifA
}

func MotifAnalysis(filename string, branches []codemlwrapper.Branch, d int, bm BlastMap) {
	_, alignedMotif, alignment, nMap := ReadAlignmentPhylip(strings.Replace(filename, ".phy", ".reconstructed.phy", -1), bm)
	f, err := os.Create(strings.Replace(filename, ".phy", ".motif.analysis.txt", -1))
	if err != nil {
		log.Panicln(err)
	}
	defer f.Close()
	writer := bufio.NewWriter(f)
	writer.WriteString("Origin\tTarget\tPosition_Source\tOrigin_Aligned_Start\tOrigin_Aligned_End\tOrigin_Seq\tTarget_Aligned_Start\tTarget_Aligned_End\tTarget_Seq\tExpanded_Origin_Start\tExpanded_Origin_End\tExpanded_Origin_Seq\tExpanded_Target_Start\tExpanded_Target_End\tExpanded_Target_seq\tMotif_Status\n")

	for _, b := range branches {
		conservedChecked := make(map[int]bool)
		if v, ok := alignedMotif[b.Origin]; ok {
			for _, e := range v {
				if _, ok:= conservedChecked[e[0]]; !ok {
					if m, ok := nMap[b.Target]; ok {
						if n, ok := m[e[0]]; ok {
							conservedChecked[e[0]] = true
							startOrigin, endOrigin := ExpandSequence(10,10, alignment.Alignment[b.Origin], e[0], e[0]+1, "-")
							startTarget, endTarget := ExpandSequence(10,10, alignment.Alignment[b.Target], e[0], e[0]+1, "-")

							writer.WriteString(b.Origin+"\t"+b.Target+"\t"+b.Origin+"\t"+strconv.Itoa(e[0])+"\t"+strconv.Itoa(e[1])+"\t"+alignment.Alignment[b.Origin][n[0]:n[1]]+"\t"+strconv.Itoa(e[0])+"\t"+strconv.Itoa(e[1])+"\t"+alignment.Alignment[b.Target][n[0]:n[1]]+"\t"+strconv.Itoa(startOrigin)+"\t"+strconv.Itoa(endOrigin)+"\t"+alignment.Alignment[b.Origin][startOrigin:endOrigin]+"\t"+strconv.Itoa(startTarget)+"\t"+strconv.Itoa(endTarget)+"\t"+alignment.Alignment[b.Target][startTarget:endTarget]+"\t"+"conserved\n")
						} else {

							start := e[0]
							if e[0] -d <= 0 {
								start = 0
							} else {
								start = e[0] -d
							}
							end := e[0] +1
							if end +d >= alignment.Length {
								end = alignment.Length
							} else {
								end = e[0] + d
							}
							s := alignment.Alignment[b.Target][start:end]
							result := nRegex.FindAllStringIndex(s, -1)
							if len(result) > 0 {
								for _, a := range result {
									if n, ok := m[a[0] + start]; ok {
										conservedChecked[e[0]] = true
										startOrigin, endOrigin := ExpandSequence(10,10, alignment.Alignment[b.Origin], e[0], e[0]+1, "-")
										startTarget, endTarget := ExpandSequence(10,10, alignment.Alignment[b.Target], n[0], n[0]+1, "-")

										writer.WriteString(b.Origin+"\t"+b.Target+"\t"+b.Origin+"\t"+strconv.Itoa(e[0])+"\t"+strconv.Itoa(e[1])+"\t"+alignment.Alignment[b.Origin][e[0]:e[1]]+"\t"+strconv.Itoa(n[0])+"\t"+strconv.Itoa(n[1])+"\t"+alignment.Alignment[b.Target][n[0]:n[1]]+"\t"+strconv.Itoa(startOrigin)+"\t"+strconv.Itoa(endOrigin)+"\t"+alignment.Alignment[b.Origin][startOrigin:endOrigin]+"\t"+strconv.Itoa(startTarget)+"\t"+strconv.Itoa(endTarget)+"\t"+alignment.Alignment[b.Target][startTarget:endTarget]+"\t"+"potentially conserved\n")
										break
									}
								}
							} else {
								startOrigin, endOrigin := ExpandSequence(10,10, alignment.Alignment[b.Origin], e[0], e[0]+1, "-")
								startTarget, endTarget := ExpandSequence(10,10, alignment.Alignment[b.Target], e[0], e[0]+1, "-")
								writer.WriteString(b.Origin+"\t"+b.Target+"\t"+b.Origin+"\t"+strconv.Itoa(e[0])+"\t"+strconv.Itoa(e[1])+"\t"+alignment.Alignment[b.Origin][e[0]:e[1]]+"\t"+strconv.Itoa(e[0])+"\t"+strconv.Itoa(e[1])+"\t"+alignment.Alignment[b.Target][e[0]:e[1]]+"\t"+strconv.Itoa(startOrigin)+"\t"+strconv.Itoa(endOrigin)+"\t"+alignment.Alignment[b.Origin][startOrigin:endOrigin]+"\t"+strconv.Itoa(startTarget)+"\t"+strconv.Itoa(endTarget)+"\t"+alignment.Alignment[b.Target][startTarget:endTarget]+"\t"+"loss\n")
							}
						}
					}
				}

			}
		}
		if v2, ok := alignedMotif[b.Target]; ok {
			for _, e := range v2 {
				if _, ok:= conservedChecked[e[0]]; !ok {
					if m, ok := nMap[b.Origin]; ok {
						if n, ok := m[e[0]]; ok {
							conservedChecked[e[0]] = true
							startOrigin, endOrigin := ExpandSequence(10,10, alignment.Alignment[b.Origin], n[0], n[0]+1, "-")
							startTarget, endTarget := ExpandSequence(10,10, alignment.Alignment[b.Target], e[0], e[0]+1, "-")
							writer.WriteString(b.Origin+"\t"+b.Target+"\t"+b.Target+"\t"+strconv.Itoa(n[0])+"\t"+strconv.Itoa(n[1])+"\t"+alignment.Alignment[b.Origin][n[0]:n[1]]+"\t"+strconv.Itoa(e[0])+"\t"+strconv.Itoa(e[1])+"\t"+alignment.Alignment[b.Target][e[0]:e[1]]+"\t"+strconv.Itoa(startOrigin)+"\t"+strconv.Itoa(endOrigin)+"\t"+alignment.Alignment[b.Origin][startOrigin:endOrigin]+"\t"+strconv.Itoa(startTarget)+"\t"+strconv.Itoa(endTarget)+"\t"+alignment.Alignment[b.Target][startTarget:endTarget]+"\t"+"conserved\n")
						} else {
							start := e[0]
							if e[0] -d <= 0 {
								start = 0
							} else {
								start = e[0] -d
							}
							end := e[0] +1
							if end +d >= alignment.Length {
								end = alignment.Length
							} else {
								end = e[1] + d
							}
							s := alignment.Alignment[b.Origin][start:end]
							result := nRegex.FindAllStringIndex(s, -1)
							if len(result) > 0 {
								for _, a := range result {
									if n, ok := m[a[0] + start]; ok {
										conservedChecked[e[0]] = true
										startOrigin, endOrigin := ExpandSequence(10,10, alignment.Alignment[b.Origin], n[0], n[0]+1, "-")
										startTarget, endTarget := ExpandSequence(10,10, alignment.Alignment[b.Target], e[0], e[0]+1, "-")
										writer.WriteString(b.Origin+"\t"+b.Target+"\t"+b.Target+"\t"+strconv.Itoa(n[0])+"\t"+strconv.Itoa(n[1])+"\t"+alignment.Alignment[b.Origin][n[0]:n[1]]+"\t"+strconv.Itoa(e[0])+"\t"+strconv.Itoa(e[1])+"\t"+alignment.Alignment[b.Target][e[0]:e[1]]+"\t"+strconv.Itoa(startOrigin)+"\t"+strconv.Itoa(endOrigin)+"\t"+alignment.Alignment[b.Origin][startOrigin:endOrigin]+"\t"+strconv.Itoa(startTarget)+"\t"+strconv.Itoa(endTarget)+"\t"+alignment.Alignment[b.Target][startTarget:endTarget]+"\t"+"potentially conserved\n")
										break
									}
								}
							} else {
								startOrigin, endOrigin := ExpandSequence(10,10, alignment.Alignment[b.Origin], e[0], e[0]+1, "-")
								startTarget, endTarget := ExpandSequence(10,10, alignment.Alignment[b.Target], e[0], e[0]+1, "-")
								writer.WriteString(b.Origin+"\t"+b.Target+"\t"+b.Target+"\t"+strconv.Itoa(e[0])+"\t"+strconv.Itoa(e[1])+"\t"+alignment.Alignment[b.Origin][e[0]:e[1]]+"\t"+strconv.Itoa(e[0])+"\t"+strconv.Itoa(e[1])+"\t"+alignment.Alignment[b.Target][e[0]:e[1]]+"\t"+strconv.Itoa(startOrigin)+"\t"+strconv.Itoa(endOrigin)+"\t"+alignment.Alignment[b.Origin][startOrigin:endOrigin]+"\t"+strconv.Itoa(startTarget)+"\t"+strconv.Itoa(endTarget)+"\t"+alignment.Alignment[b.Target][startTarget:endTarget]+"\t"+"gain\n")
							}
						}
					}
				}
				//log.Println(alignment.Alignment[b.Target])
				//log.Println(alignment.Alignment[b.Target][e2[0]:e2[1]])
			}
		}
	}
	writer.Flush()
}

func CombineTree(a codemlwrapper.CodeMLCommandline) {
	p := pythoncmdwrapper.ProcessTreeCommandline{}
	p.Command = `C:\Program Files\Anaconda3\python.exe`
	p.Program = `D:\GoProject\ancestral\src\github.com\noatgnu\ancestral\processtree.py`
	p.CurrentTree = a.TreeFile
	p.Reconstructed = strings.Replace(a.TreeFile, "_tree", "_reconstructed_tree", -1)
	p.Out = strings.Replace(p.Reconstructed, ".txt", ".reconstructed.tree.txt", -1)
	err := p.Execute()
	if err != nil {
		log.Panicln(err)
	}
}

func ProcessAlignment(b BlastMap, asr bool) {
	log.Printf("Started: Phylogeny Construction (%v)", b.FileName)
	alignmentFile := strings.Replace(b.FileName, ".filtered.fasta", ".phy", -1)
	CreateAlignment(b.FileName, alignmentFile)
	if asr == true {
		CreateTreeML(alignmentFile, 1000, "fast")
		ReformatBootstrappedTree(alignmentFile)

		ASR(alignmentFile, alignmentFile+"_phyml_tree_reformated.txt", strings.Replace(alignmentFile, ".phy", ".asr", -1), b)
	} else {
		ReadAlignmentPhylip(alignmentFile, b)
	}
}

func ReformatBootstrappedTree(alignmentFile string) {

	f, err := os.Open(alignmentFile + "_phyml_tree.txt")
	if err != nil {
		log.Panicln(err)
	}
	reader := bufio.NewReader(f)
	w, err := os.Create(alignmentFile + "_phyml_tree_reformated.txt")
	if err != nil {
		log.Panicln(err)
	}
	writer := bufio.NewWriter(w)
	r, _, err := reader.ReadLine()
	if err != nil {
		log.Panicln(err)
	}
	line := strings.TrimSpace(string(r[:]))
	log.Printf("Reformatting Bootstrapped Tree %v", line)
	matches := treeRegex.FindAllStringSubmatchIndex(line, -1)
	replaced := make(map[string]bool)
	if len(matches) > 0 {
		for _, m := range matches {
			r := line[m[0]:m[1]]
			if _, ok := replaced[r]; !ok {
				replaced[r] = true
			}
		}
	}
	for k := range replaced {
		line = strings.Replace(line, k, "):", -1)
	}
	writer.WriteString(line + "\n")
	writer.Flush()
	f.Close()
	w.Close()
}

func ExpandSequence(numberForward int, numberBackward int, sequence string, start int, end int, ignoreCharacter string) (startPosition int, endPosition int) {
	forwardOffset := 0
	backwardOffset := 0
	lSeq := len(sequence)
	for {
		if numberForward == 0 && numberBackward == 0 {
			break
		}
		if numberForward > 0 {
			l := end + forwardOffset
			if l >= lSeq {
				numberForward = 0
			} else {
				if sequence[l:l+1] != ignoreCharacter {
					numberForward --
				}
				forwardOffset ++
			}

		}
 		if numberBackward > 0 {
			backwardOffset ++
			l := start - backwardOffset
			if l <= 0 {
				numberBackward = 0
			} else {
				if sequence[l:l+1] != ignoreCharacter {
					numberBackward --
				}
			}
		}
	}
	return start-backwardOffset, end+forwardOffset
}

func RemapTopDom(seqWithGap string, topDom []blastwrapper.TopDom, filename string) (reMapped []blastwrapper.TopDom) {
	log.Println("Started: Topological Domain Remap")
	o, err := os.Create(filename)
	if err != nil {
		log.Panicln(err)
	}
	log.Println(seqWithGap)
	writer := csv.NewWriter(o)
	writer.Comma = '\t'
	writer.Write([]string{"Origin_Start", "Origin_End", "Remapped_Start", "Remapped_End", "Type"})
	g := gappedRegex.FindAllStringSubmatchIndex(seqWithGap, -1)

	for _, t := range topDom {
		log.Println(t)
		gap := 0
		startCheck := false
		endCheck := false
		td := blastwrapper.TopDom{}
		sequence := seqWithGap[0:t.Start]
		log.Println(sequence)
		if strings.Contains(sequence, "-") {
			for _, v := range g {
				seq := seqWithGap[v[0]:v[1]]
				log.Println(seq)
				gap += strings.Count(seq, "-")
				if startCheck == false {
					if (v[0] <= (t.Start+gap-1)) && (v[1] >(t.Start+gap)) {
						td.Start = t.Start+gap
						startCheck = true
					}
				}
				if endCheck == false {
					if (v[0] <= (t.Stop+gap)) && (v[1] >= (t.Stop+gap)) {
						td.Stop = t.Stop+gap
						endCheck = true
					}
				}
				if (startCheck == true) && (endCheck == true) {
					break
				}
			}
		} else {
			td.Start = t.Start
			if strings.Contains(seqWithGap[t.Start-1:t.Stop], "-") == true {
				for _, v := range g {
					seq := seqWithGap[v[0]:v[1]]
					gap += strings.Count(seq, "-")
					if endCheck == false {
						if (v[0] <= (t.Stop+gap)) && (v[1] >= (t.Stop+gap)) {
							td.Stop = t.Stop+gap
							endCheck = true
							break
						}
					}
				}
			} else {
				td.Stop = t.Stop
			}
		}
		td.Type = t.Type
		log.Println(t, td)
		writer.Write([]string{strconv.Itoa(t.Start), strconv.Itoa(t.Stop), strconv.Itoa(td.Start), strconv.Itoa(td.Stop), t.Type})
		reMapped = append(reMapped, td)
	}
	writer.Flush()
	defer log.Println("Finished: Topological Domain Remap")
	defer o.Close()

	return reMapped
}