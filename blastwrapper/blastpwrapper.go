package blastwrapper

import (
	"errors"
	"os/exec"
	"strconv"
	"log"
)

type NcbiBlastpCommandline struct {
	Command    string
	Query      string
	QueryLoc   string
	DB         string
	EValue     float64
	NumThreads int
	Out        string
	OutFmt     string
	CommandLine
}

func (b *NcbiBlastpCommandline) Execute() (err error) {
	commandArray, err := b.CommandBuild()
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

func (b *NcbiBlastpCommandline) CommandBuild() (commandArray []string, err error) {
	if b.Command != "" {
		commandArray = append(commandArray, b.Command)
	} else {
		commandArray = append(commandArray, "blastp")
	}
	if b.Query != "" {
		commandArray = append(commandArray, "-query", b.Query)
		if b.QueryLoc != "" {
			commandArray = append(commandArray, "-query_loc", b.QueryLoc)
		}
	} else {
		return nil, errors.New("text: need query input")
	}

	if b.DB != "" {
		commandArray = append(commandArray, "-db", b.DB)
	} else {
		return nil, errors.New("text: need db name")
	}

	if b.EValue != 0 {
		commandArray = append(commandArray, "-evalue", strconv.FormatFloat(b.EValue, 'f', -1, 64))
	} else {
		commandArray = append(commandArray, "-evalue", "0.001")
	}

	if b.NumThreads != 0 {
		commandArray = append(commandArray, "-num_threads", strconv.Itoa(b.NumThreads))
	}

	if b.Out != "" {
		commandArray = append(commandArray, "-out", b.Out)
	}

	if b.OutFmt != "" {
		commandArray = append(commandArray, "-outfmt", b.OutFmt)
	} else {
		commandArray = append(commandArray, "-outfmt", `"6 qseqid sgi sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"`)
	}
	return commandArray, nil
}
