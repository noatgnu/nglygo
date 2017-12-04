package clustalowrapper

import (
	"github.com/noatgnu/ancestral/blastwrapper"
	"errors"
	"os/exec"
	"log"
)

type ClustalOCommandline struct {
	Command string
	In string
	Out string
	OutFmt string
	blastwrapper.CommandLine
}

func (c *ClustalOCommandline) Execute() (err error) {
	commandArray, err := c.CommandBuild()
	if err != nil {
		return err
	}
	if commandArray != nil {
		cmd := exec.Command(commandArray[0], commandArray[1:]...)
		log.Println(cmd.Args)
		_ = cmd.Run()
	}
	return err
}

func (c *ClustalOCommandline) CommandBuild() (commandArray []string, err error) {
	if c.Command != "" {
		commandArray = append(commandArray, c.Command)
	} else {
		commandArray = append(commandArray, `clustalo`)
	}
	if c.In != "" {
		commandArray = append(commandArray, "-i", c.In)
	} else {
		return nil, errors.New("text: Require Input")
	}
	if c.Out != "" {
		commandArray = append(commandArray, "-o", c.Out)
	} else {
		return nil, errors.New("text: Require Output")
	}
	if c.OutFmt != "" {
		commandArray = append(commandArray, c.OutFmt)
	}

	return commandArray, nil
}
