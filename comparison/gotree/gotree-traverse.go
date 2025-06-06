package main

import (
	"fmt"
	"os"
	"time"

	"github.com/evolbioinfo/gotree/io/newick"
	"github.com/evolbioinfo/gotree/tree"
)

func main() {
	var t *tree.Tree
	var err error
	var f *os.File
	if f, err = os.Open(os.Args[1]); err != nil {
		panic(err)
	}
	var start = time.Now()
	t, err = newick.NewParser(f).Parse()
	if err != nil {
		panic(err)
	}
	t.PostOrder(func(cur *tree.Node, prev *tree.Node, e *tree.Edge) bool {
		if e != nil {
			// fmt.Print("%v", e.Length()))
		}
		if cur.Tip() {
			fmt.Println(cur.Name())
		} else {
		}
		return true
	})

	elapsed := time.Since(start).Seconds()
	fmt.Printf("\nInternal time: %v\n", elapsed)

}
