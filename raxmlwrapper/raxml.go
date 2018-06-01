package raxmlwrapper

import (
	"github.com/noatgnu/ancestral/blastwrapper"
	"os/exec"
	"log"
	"strconv"
	"errors"
	"path/filepath"
)


type RaxMLCommandline struct {
	Command        string
	In             string
	Model string
	Bootstrap int
	Out string
	blastwrapper.CommandLine
}



func (p *RaxMLCommandline) Execute() (err error) {
	commandArray, err := p.CommandBuild()
	if err !=nil {
		return err
	}
	if commandArray !=nil {
		print(commandArray)
		cmd := exec.Command(commandArray[0], commandArray[1:]...)
		log.Printf("Started: Creating Phylogenetic Tree (%v)", p.In)
		cmd.Run()
		log.Printf("Finished: Creating Phylogenetic Tree (%v)", p.In)
	}
	return err
}

func (p *RaxMLCommandline) CommandBuild() (commandArray []string, err error) {
	if p.Command != "" {
		commandArray = append(commandArray, p.Command)
	} else {
		commandArray = append(commandArray, "raxML")
	}
	if p.In != "" {
		commandArray = append(commandArray, "-s", p.In)
	} else {
		return nil, errors.New("require Input")
	}

	if p.Bootstrap != 0 {
		commandArray = append(commandArray, "--bootstop-perms="+strconv.Itoa(p.Bootstrap))
	} else {
		commandArray = append(commandArray, "--bootstop-perms=100")
	}

	if p.Model != "" {
		commandArray = append(commandArray, "-m", p.Model)
	} else {
		commandArray = append(commandArray, "-m","PROTCATWAG")
	}

	commandArray = append(commandArray, "-slow")
	if p.Out != "" {
		folder, file := filepath.Split(p.Out)
		commandArray = append(commandArray,"-w",folder,"-n", file)
	} else {
		return nil, errors.New("require Output")
	}


	return commandArray, nil
}

