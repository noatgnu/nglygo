package blastwrapper

import (
	"errors"
	"log"
	"os/exec"
)

type MakeBlastDBCommandline struct {
	Command string
	In      string
	DB      string
	Out     string
	CommandLine
}

func (m *MakeBlastDBCommandline) Execute() (err error) {
	commandArray, err := m.CommandBuild()
	if err != nil {
		return err
	}
	if commandArray != nil {
		cmd := exec.Command(commandArray[0], commandArray[1:]...)
		log.Println(cmd.Args)
		err := cmd.Run()
		if err != nil {
			return err
		}
	} else {
		return errors.New("text: need parameters")
	}
	return err
}

func (m *MakeBlastDBCommandline) CommandBuild() (commandArray []string, err error) {
	if m.Command != "" {
		commandArray = append(commandArray, m.Command)
	} else {
		commandArray = append(commandArray, "makeblastdb")
	}
	if m.In != "" {
		commandArray = append(commandArray, "-in", m.In)
	} else {
		return nil, errors.New("text: need input file")
	}
	if m.DB != "" {
		commandArray = append(commandArray, "-dbtype", m.DB)
	} else {
		return nil, errors.New("text: need db type")
	}
	if m.Out != "" {
		commandArray = append(commandArray, "-out", m.Out)
	} else {
		return nil, errors.New("text: need output file")
	}
	commandArray = append(commandArray, "-parse_seqids")

	return commandArray, nil
}
