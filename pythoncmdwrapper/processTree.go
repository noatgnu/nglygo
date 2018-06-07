package pythoncmdwrapper

import (
	"github.com/noatgnu/ancestral/blastwrapper"
	"errors"
	"os/exec"
	"log"
	"bytes"
)

type ProcessTreeCommandline struct {
	Command string
	Program string
	CurrentTree string
	Reconstructed string
	Out string
	blastwrapper.CommandLine
}

func (p *ProcessTreeCommandline) Execute() (err error) {
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

func (p *ProcessTreeCommandline) CommandBuild() (commandArray []string, err error) {
	if p.Command != "" {
		commandArray = append(commandArray, p.Command)
	} else {
		commandArray = append(commandArray, "python")
	}
	if p.Program != "" {
		commandArray = append(commandArray, p.Program)
	} else {
		commandArray = append(commandArray, "processTree.py")
	}
	if p.CurrentTree != "" {
		commandArray = append(commandArray, "-ct", p.CurrentTree)
	} else {
		return nil, errors.New("require current tree")
	}
	if p.Reconstructed != "" {
		commandArray = append(commandArray, "-asrt", p.Reconstructed)
	} else {
		return nil, errors.New("require reconstructed tree")
	}
	if p.Out != "" {
		commandArray = append(commandArray, "-o", p.Out)
	} else {
		return nil, errors.New("require output filename")
	}
	return commandArray, nil
}
