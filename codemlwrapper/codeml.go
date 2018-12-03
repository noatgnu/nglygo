package codemlwrapper

import (
	"github.com/noatgnu/ancestral/blastwrapper"
	"errors"
	"os"
	"bufio"
	"fmt"
	"os/exec"
	"log"
	"path/filepath"
	"bytes"
	"io"
	"strings"
	"regexp"
	"time"
)

type Branch struct{
	Origin string
	Target string
	Changes [][]string
}

type BranchAnalysis struct{
	Results [][]string `json:"results"`
}

const branchRegexPat = `(\d+)..(\d+)`
var branchRegex = regexp.MustCompile(branchRegexPat)

const changeRegexPat = `([0-9]*)\s(.)\s([0-9]*\.?[0-9]*)\s->\s(.)`
var changeRegex = regexp.MustCompile(changeRegexPat)

type CodeMLCommandline struct {
	Command string
	Ctl string
	SeqFile string
	OutFile string
	TreeFile string
	Noisy int
	Verbose int
	Runmode int
	SeqType int
	Clock int
	AADist int
	AARateFile string
	Model int
	ICode int
	Mgene int
	FixAlpha int
	Alpha string
	MAlpha int
	NCatG int
	GetSE int
	RateAncestor int
	SmallDiff string
	CleanData int
	Method int
	NData int
	blastwrapper.CommandLine
}

func (b *Branch) ToFile(writer *bufio.Writer) {
	for _, v := range b.Changes {
		writer.WriteString(fmt.Sprintf("%v\t%v\t%v\t%v\t%v\t%v\n", b.Origin, b.Target, v[0], v[1], v[2], v[3]))
	}
}

func (c *CodeMLCommandline) Execute() (err error) {
	commandArray, err := c.CommandBuild()
	if err != nil {
		return err
	}
	if commandArray != nil {
		log.Println(commandArray)
		cmd := exec.Command(commandArray[0], commandArray[1:]...)
		dir, _ := filepath.Split(c.SeqFile)
		var stderr bytes.Buffer
		cmd.Stderr = &stderr
		cmd.Dir = dir
		log.Printf("Started: Ancestral Sequence Reconstruction (%v)", cmd.Dir)
		err = cmd.Run()
		if err != nil {
			log.Panicln(c.SeqFile + stderr.String())
		}
		log.Printf("Finished: Ancestral Sequence Reconstruction (%v)", cmd.Dir)
	}
	return err
}

func (c *CodeMLCommandline) CommandBuild() (commandArray []string, err error) {
	if c.Command != "" {
		commandArray = append(commandArray, c.Command)
	} else {
		commandArray = append(commandArray, "codeml")
	}
	if c.Ctl != "" {
		commandArray = append(commandArray, c.Ctl)
	} else {
		return nil, errors.New("required configuration")
	}
	return commandArray, nil
}

func (c *CodeMLCommandline) BuildCtl(fileName string) (err error) {
	f, err := os.Create(fileName)
	writer := bufio.NewWriter(f)
	if c.SeqFile != "" {
		writer.WriteString(fmt.Sprintf("seqfile = %v\n", c.SeqFile))
	} else {
		return errors.New("require input sequence")
	}
	if c.OutFile != "" {
		writer.WriteString(fmt.Sprintf("outfile = %v\n", c.OutFile))
	} else {
		return errors.New("require output file")
	}
	if c.TreeFile != "" {
		writer.WriteString(fmt.Sprintf("treefile = %v\n", c.TreeFile))
	} else {
		return errors.New("require tree file")
	}
	if c.AARateFile != "" {
		writer.WriteString(fmt.Sprintf("aaRatefile = %v\n", c.AARateFile))
	} else {
		return errors.New("require AA substitution rate file")
	}
	writer.WriteString(
		fmt.Sprintf("noisy = %v\nverbose = %v\nrunmode = %v\nseqtype = %v\nclock = %v\naaDist = %v\nmodel = %v\nicode = %v\nMgene = %v\nfix_alpha = %v\nMalpha = %v\nncatG = %v\ngetSE = %v\nRateAncestor = %v\ncleandata = %v\nmethod = %v\nndata = %v\n", c.Noisy, c.Verbose, c.Runmode, c.SeqType, c.Clock, c.AADist, c.Model, c.ICode, c.Mgene, c.FixAlpha, c.MAlpha, c.NCatG, c.GetSE, c.RateAncestor, c.CleanData, c.Method, c.NData))
	if c.SmallDiff != "" {
		writer.WriteString(fmt.Sprintf("Small_Diff = %v\n", c.SmallDiff))
	} else {
		return errors.New("require smalldiff parameter")
	}
	if c.Alpha != "" {
		writer.WriteString(fmt.Sprintf("alpha = %v", c.Alpha))
	} else {
		return errors.New("require alpha parameter")
	}

	c.Ctl = fileName
	defer writer.Flush()
	return err
}

// Reading supplementary result from CodeML ancestral sequence reconstruction operation. Parse out tree and branches mutation
// information. Create a _branch file to store the result. If no reconstructed tree could be found, the process will pause
// for 2 seconds and the file would be re read. Return all collected branches from the reconstructed information.
func (c *CodeMLCommandline) ReadSupplemental(path string) (branches []Branch) {

	f, err := os.Open(path)
	if err != nil {
		log.Panicln(err)
	}
	log.Printf("Reading rst file %v", path)
	buff := bufio.NewReader(f)
	bf, err := os.Create(c.SeqFile+"_branch")
	if err != nil {
		log.Panicln(err)
	}

	writer := bufio.NewWriter(bf)
	writer.WriteString("Origin\tTarget\tPosition\tOrigin_Residue\tStats\tTarget_Residue\n")
	//branches = make(map[string]*Branch)
	lineChan := make(chan string)
	go func() {
		for {
			r, err := buff.ReadString('\n')
			if err == io.EOF {
				break
			}
			if err != nil {
				log.Panicln(err)
			}
			s := strings.TrimSpace(string(r))
			s = strings.TrimSpace(s)
			lineChan <- s
		}
		close(lineChan)
	}()
	var tree string
	var RA bool
	for s := range lineChan {
		if s != "" {
			if strings.HasPrefix(s,"tree with node labels") {
				treeFileName := strings.Replace(c.TreeFile, "_tree", "_reconstructed_tree", -1)
				f, err := os.Create(treeFileName)
				if err != nil {
					log.Panicln(err)
				}
				writer := bufio.NewWriter(f)
				tree = ReadNextLine(lineChan)
				writer.WriteString(tree+"\n")
				writer.Flush()
				f.Close()
				log.Printf("Wrote Reconstructed Tree %v", treeFileName)

			} else if strings.HasPrefix(s, "List of extant and reconstructed sequences") {
				log.Printf("Getting Reconstructed Alignment %v", c.SeqFile)
				RA = true
				WriteReconstructedAlignment(c.SeqFile ,lineChan)

			} else if strings.HasPrefix(s, "Branch") {
				var branch Branch
				result := branchRegex.FindAllStringSubmatch(s, 2)
				branch.Origin = result[0][1]
				branch.Target = result[0][2]
				ReadBranch(lineChan, &branch, writer)
				branches = append(branches, branch)
			}
		}
	}

	writer.Flush()
	f.Close()
	bf.Close()
	if tree == "" || RA != true {
		time.Sleep(2)
		c.ReadSupplemental(path)
	}

	return branches
}

func ReadNextLine(lineChan chan string) string {
	for s := range lineChan {
		return s
		break
	}
	return ""

}

// Write out reconstructed alignment from supplementary file.
func WriteReconstructedAlignment(filename string, lineChan chan string) {
	af, err := os.Create(strings.Replace(filename, ".phy", ".reconstructed.phy", -1))
	if err != nil {
		log.Panicln(err)
	}

	writer := bufio.NewWriter(af)
	alignmentStarted := false
	started := false
	for s := range lineChan {
		if s != ""  {
			if started && !alignmentStarted {
				writer.WriteString(" "+s+"\n")
			} else {
				writer.WriteString(s[0:18]+strings.Replace(s[18:], " ", "", -1)+"\n")
			}
		} else {
			if !started {
				started = true
			} else {
				if !alignmentStarted {
					alignmentStarted = true
				}  else {
					break
				}
			}

		}
	}
	writer.Flush()
	af.Close()
	log.Printf("Wrote Reconstructed Alignment for %v", filename)
}

// Read string from channel and identify branches information and write out to output writer.
func ReadBranch(lineChan chan string, branch *Branch, writer *bufio.Writer) {
	enterBranch := false
	emptyLine := 0
	for s := range lineChan {
		if s != "" {
			enterBranch = true
			changes := changeRegex.FindAllStringSubmatch(s, -1)
			if changes != nil {
				writer.WriteString(fmt.Sprintf("%v\t%v\t%v\t%v\t%v\t%v\n", branch.Origin, branch.Target, changes[0][1], changes[0][2], changes[0][3], changes[0][4]))
				branch.Changes = append(branch.Changes, changes[0][1:])
			}
		} else if s == "" {
			emptyLine ++
			if enterBranch {
				break
			}
			if emptyLine == 3 {
				break
			}
		}
	}

}