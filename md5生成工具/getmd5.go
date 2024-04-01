package main
import (
    "crypto/md5"
    "fmt"
    "io"
    _"io/ioutil"
    "os"
    "path/filepath"
)
func md5sum(filepath string) (string, error) {
    file, err := os.Open(filepath)
    if err != nil {
        return "", err
    }
    defer file.Close()
    hasher := md5.New()
    if _, err := io.Copy(hasher, file); err != nil {
        return "", err
    }
    return fmt.Sprintf("%x", hasher.Sum(nil)), nil
}
func generate_md5sum_file(dirpath string) error {
    file_count := 0
    filepath.Walk(dirpath, func(path string, info os.FileInfo, err error) error {
        if err != nil {
            return err
        }
        if info.IsDir() {
            return nil
        }
        file_count++
        return nil
    })
    fp, err := os.Create(filepath.Join(dirpath, "MD5sum.txt"))
    if err != nil {
        return err
    }
    defer fp.Close()
    processed_files := 0
    filepath.Walk(dirpath, func(path string, info os.FileInfo, err error) error {
        if err != nil {
            return err
        }
        if info.IsDir() {
            return nil
        }
        md5, err := md5sum(path)
        if err != nil {
            return err
        }
        relpath, err := filepath.Rel(dirpath, path)
        if err != nil {
            return err
        }
        fmt.Fprintf(fp, "%s  %s\n", md5, relpath)
        processed_files++
        progress := processed_files * 100 / file_count
        fmt.Printf("\rProcessed %d/%d files [%d%%] ", processed_files, file_count, progress)
        for i := 0; i < progress; i++ {
            fmt.Print("=")
        }
        for i := progress; i < 100; i++ {
            fmt.Print(" ")
        }
        return nil
    })
    return nil
}
func main() {
    if len(os.Args) != 3 || os.Args[1] != "--dir" {
        fmt.Fprintf(os.Stderr, "Usage: %s --dir directory\n", os.Args[0])
        os.Exit(1)
    }
    dirpath := os.Args[2]
    err := generate_md5sum_file(dirpath)
    if err != nil {
        fmt.Fprintf(os.Stderr, "Error: %v\n", err)
        os.Exit(1)
    }
    fmt.Println("\nMD5sum generation completed successfully!")
}
