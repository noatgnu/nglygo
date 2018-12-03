package pythoncmdwrapper

import (
	"bytes"
	"errors"
	"github.com/noatgnu/ancestral/blastwrapper"
	"log"
	"os/exec"
)

type PruneTreeCommandline struct {
	Command string
	Program string
	CurrentTree string
	Organism string
	Replace string
	Out string
	blastwrapper.CommandLine
}

func (p *PruneTreeCommandline) Execute() (err error) {
	commandArray, err := p.CommandBuild()
	if err != nil {
		return err
	}
	if commandArray != nil {
		log.Println(commandArray)
		cmd := exec.Command(commandArray[0], commandArray[1:]...)
		var stderr bytes.Buffer
		cmd.Stderr = &stderr
		err = cmd.Run()

		if err != nil {
			return errors.New(stderr.String())
		}
	}
	return err
}

func (p *PruneTreeCommandline) CommandBuild() (commandArray []string, err error) {
	if p.Command != "" {
		commandArray = append(commandArray, p.Command)
	} else {
		commandArray = append(commandArray, "python")
	}
	if p.Program != "" {
		commandArray = append(commandArray, p.Program)
	} else {
		commandArray = append(commandArray, "pruneTree.py")
	}
	if p.CurrentTree != "" {
		commandArray = append(commandArray, `-in`, p.CurrentTree)
	} else {
		return nil, errors.New("require current tree")
	}
	if p.Organism != "" {
		commandArray = append(commandArray, `-org`, p.Organism)
	} else {
		return nil, errors.New("require organism to be retained")
	}
	if p.Replace != "" {
		commandArray = append(commandArray, `-r`, p.Replace)
	}

	if p.Out != "" {
		commandArray = append(commandArray, `-o`, p.Out)
	} else {
		return nil, errors.New("require output filename")
	}
	return commandArray, nil
}
