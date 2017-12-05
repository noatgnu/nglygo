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
)

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
	blastwrapper.CommandLine
}

func (c *CodeMLCommandline) Execute() (err error) {
	commandArray, err := c.CommandBuild()
	if err != nil {
		return err
	}
	if commandArray != nil {
		cmd := exec.Command(commandArray[0], commandArray[1:]...)
		dir, _ := filepath.Split(c.SeqFile)
		var stderr bytes.Buffer
		cmd.Stderr = &stderr
		cmd.Dir = dir
		log.Printf("Started: Ancestral Sequence Reconstruction (%v)", c.SeqFile)
		err = cmd.Run()
		if err != nil {
			log.Panicln(stderr.String())
		}
		log.Printf("Finished: Ancestral Sequence Reconstruction (%v)", c.SeqFile)
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
		fmt.Sprintf("noisy = %v\nverbose = %v\nrunmode = %v\nseqtype = %v\nclock = %v\naaDist = %v\nmodel = %v\nicode = %v\nMgene = %v\nfix_alpha = %v\nMalpha = %v\nncatG = %v\ngetSE = %v\nRateAncestor = %v\ncleandata = %v\nmethod = %v\n", c.Noisy, c.Verbose, c.Runmode, c.SeqType, c.Clock, c.AADist, c.Model, c.ICode, c.Mgene, c.FixAlpha, c.MAlpha, c.NCatG, c.GetSE, c.RateAncestor, c.CleanData, c.Method))
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

func (c *CodeMLCommandline) ReadSupplemental() {
	dir, _ := filepath.Split(c.SeqFile)
	f, err := os.Open(filepath.Join(dir, "rst"))
	if err != nil {
		log.Panicln(err)
	}
	buff := bufio.NewReader(f)
	currentTree := False
	for {
		r, err := buff.ReadString('\n')
		if err != nil {
			if err == io.EOF {
				break
			}
			log.Panicln(err)
		}
		s := strings.TrimSpace(r)
		if s == "Ancestral reconstruction by AAML." {

		}
	}
}
