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
	var nbtips int
	if f, err = os.Open(os.Args[1]); err != nil {
		panic(err)
	}
	t, err = newick.NewParser(f).Parse()
	nbtips = len(t.Tips())
	if err != nil {
		panic(err)
	}
	var start = time.Now()

	t, err = tree.RandomYuleBinaryTree(nbtips, true)
	//t, err = tree.RandomBalancedBinaryTree(depth, rooted)
	//t, err = tree.RandomUniformBinaryTree(nbtips, rooted)
	//t, err = tree.RandomCaterpilarBinaryTree(nbtips, rooted)
	//t, err = tree.StarTree(nbtips)

	if err != nil {
		panic(err)
	}
	fmt.Println(t.Newick())

	elapsed := time.Since(start).Seconds()
	fmt.Printf("\nInternal time: %v\n", elapsed)

}
