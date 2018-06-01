package fasttreecmdwrapper

import (
	"github.com/noatgnu/ancestral/blastwrapper"
	"os/exec"
	"log"
	"errors"
	"strconv"
)

type FastMLCommandline struct {
	Command        string
	In             string
	Model string
	Bootstrap int
	Out string
	blastwrapper.CommandLine
}



func (p *FastMLCommandline) Execute() (err error) {
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

func (p *FastMLCommandline) CommandBuild() (commandArray []string, err error) {
	if p.Command != "" {
		commandArray = append(commandArray, p.Command)
	} else {
		commandArray = append(commandArray, "FastTree")
	}

	if p.Bootstrap != 0 {
		commandArray = append(commandArray, "-boot", strconv.Itoa(p.Bootstrap))
	}

	if p.Model != "" {
		commandArray = append(commandArray, p.Model)
	} else {
		commandArray = append(commandArray, "-wag")
	}

	commandArray = append(commandArray, "-slow")
	if p.Out != "" {
		commandArray = append(commandArray, "-out", p.Out)
	} else {
		return nil, errors.New("require Output")
	}
	if p.In != "" {
		commandArray = append(commandArray, p.In)
	} else {
		return nil, errors.New("require Input")
	}

	return commandArray, nil
}
