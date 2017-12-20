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
)

const nRegexPat = `N`
var nRegex = regexp.MustCompile(nRegexPat)

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

func ReadAlignmentPhylip(filename string) (map[string][][]int, map[string][][]int, alignio.Alignment, map[string]map[int][]int) {
	a := alignio.ReadPhylip(filename)
	f, err := os.Create(strings.Replace(filename, ".phy", ".motifs.txt", -1))
	if err != nil {
		log.Panicln(err)
	}
	defer f.Close()
	writer := bufio.NewWriter(f)
	writer.WriteString("Entry\tStart_Position\tEnd_Position\tMotif\tStart_Position_With_Gap\tEnd_Position_With_Gap\tMotif_With_Gap\n")

	motifMap := make(map[string][][]int)
	alignedMotifMap := make(map[string][][]int)
	nResidueAlignedMotifMap := make(map[string]map[int][]int)
	for k, v := range a.Alignment {
		nResidueAlignedMotifMap[k] = make(map[int][]int)
		gap := nglycan.MotifParseStringWithGap(v, -1)
		noGapSeq := strings.Replace(v, "-", "", -1)
		noGap := nglycan.MotifParseStringWithGap(noGapSeq, -1)
		motifMap[k] = noGap

		alignedMotifMap[k] = gap
		for i := range gap {
			nResidueAlignedMotifMap[k][gap[i][0]] = gap[i]
			writer.WriteString(fmt.Sprintf("%v\t%v\t%v\t%v\t%v\t%v\t%v\n", k, noGap[i][0]+1, noGap[i][1], noGapSeq[noGap[i][0]:noGap[i][1]], gap[i][0]+1, gap[i][1], v[gap[i][0]:gap[i][1]]))
		}
	}
	writer.Flush()
	return motifMap, alignedMotifMap, a, nResidueAlignedMotifMap
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
	a.Verbose = 2
	a.NData = 1
	a.BuildCtl(strings.Replace(seq, ".phy", ".ctl", -1))
	err := a.Execute()
	if err != nil {
		log.Panicln(err)
	}
	branches := a.ReadSupplemental()
	MotifAnalysis(seq, branches, 0)
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

func MotifAnalysis(filename string, branches []codemlwrapper.Branch, d int) {
	_, alignedMotif, alignment, nMap := ReadAlignmentPhylip(strings.Replace(filename, ".phy", ".reconstructed.phy", -1))
	f, err := os.Create(strings.Replace(filename, ".phy", ".motif.analysis.txt", -1))
	if err != nil {
		log.Panicln(err)
	}
	defer f.Close()
	writer := bufio.NewWriter(f)
	writer.WriteString("Origin\tTarget\tPosition_Source\tOrigin_Aligned_Start\tOrigin_Aligned_End\tOrigin_Seq\tTarget_Aligned_Start\tTarget_Aligned_End\tTarget_Seq\tMotif_Status\tExpanded_Origin_Start\tExpanded_Origin_End\tExpanded_Origin_Seq\tExpanded_Target_Start\tExpanded_Target_End\tExpanded_Target_seq\n")

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
	p.Program = `C:\Users\localadmin\GoglandProjects\ancestral\src\github.com\noatgnu\ancestral\processtree.py`
	p.CurrentTree = a.TreeFile
	p.Reconstructed = strings.Replace(a.TreeFile, "_tree", "_reconstructed_tree", -1)
	p.Out = strings.Replace(p.Reconstructed, ".txt", ".reconstructed.tree.txt", -1)
	err := p.Execute()
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