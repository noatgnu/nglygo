package phymlwrapper

import (
	"github.com/noatgnu/ancestral/blastwrapper"
	"errors"
	"strconv"
)

type PhyMLCommandline struct {
	Command        string
	In             string
	DType             string
	Model string
	Alpha string
	Bootstrap int
	blastwrapper.CommandLine
}

func (p *PhyMLCommandline) Execute() (err error) {
	return err
}

func (p *PhyMLCommandline) CommandBuild() (commandArray []string, err error) {
	if p.Command != "" {
		commandArray = append(commandArray, p.Command)
	} else {
		commandArray = append(commandArray, "phyml")
	}
	if p.In != "" {
		commandArray = append(commandArray, "-i", p.In)
	} else {
		return nil, errors.New("Require Input")
	}
	if p.DType != "" {
		commandArray = append(commandArray, "-d", p.DType)
	} else {
		commandArray = append(commandArray, "-d", "aa")
	}
	if p.Model != "" {
		commandArray = append(commandArray, "-m", p.Model)
	} else {
		commandArray = append(commandArray, "-m", "WAG")
	}
	if p.Alpha != "" {
		commandArray = append(commandArray, "-a", p.Alpha)
	} else {
		commandArray = append(commandArray, "-a", "e")
	}
	if p.Bootstrap != 0 {
		commandArray = append(commandArray, "-b", strconv.Itoa(p.Bootstrap))
	} else {
		commandArray = append(commandArray, "-b", "0")
	}
	return commandArray, nil
}