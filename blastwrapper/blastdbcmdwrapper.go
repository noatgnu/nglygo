package blastwrapper

import (
	"errors"
	"log"
	"os/exec"
	"bytes"
)

type BlastDBCMDCommandline struct {
	Command        string
	In             string
	DB             string
	DBType         string
	Out            string
	FilterOrganism []string
	PartialCheck   bool
	CommandLine
}

func (b *BlastDBCMDCommandline) Execute() (err error) {
	commandArray, err := b.CommandBuild()
	if err != nil {
		return err
	}

	if commandArray != nil {
		cmd := exec.Command(commandArray[0], commandArray[1:]...)
		var stderr bytes.Buffer
		cmd.Stderr = &stderr
		log.Println(cmd.Args)
		cmd.Run()

	} else {
		return errors.New("text: need parameters")
	}
	return err
}

func (b *BlastDBCMDCommandline) CommandBuild() (commandArray []string, err error) {
	if b.Command != "" {
		commandArray = append(commandArray, b.Command)
	} else {
		commandArray = append(commandArray, "blastdbcmd")
	}

	if b.DB != "" {
		commandArray = append(commandArray, "-db", b.DB)
	} else {
		return nil, errors.New("text: need db type")
	}
	if b.DBType != "" {
		commandArray = append(commandArray, "-dbtype", b.DBType)
	}
	commandArray = append(commandArray, "-entry_batch")
	if b.In != "" {
		commandArray = append(commandArray, b.In)
	} else {
		return nil, errors.New("text: need input file")
	}
	if b.Out != "" {
		commandArray = append(commandArray, "-out", b.Out)
	} else {
		return nil, errors.New("text: need output file")
	}
	return commandArray, nil
}

