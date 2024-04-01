import os
import hashlib
import tqdm
import pandas as pd
from concurrent.futures import ProcessPoolExecutor
import argparse
from celescope.tools.utils import add_log

def generate_md5_for_files(directory):
    file_md5_records = []

    for root, dirs, files in os.walk(directory):
        for file in tqdm.tqdm(files, ncols = 50):
                        
            file_path = os.path.join(root, file)
            md5_hash = hashlib.md5()

            try :
                with open(file_path, 'rb') as f:
                    for chunk in iter(lambda: f.read(8192), b''):
                        md5_hash.update(chunk)
            except PermissionError:
                print(f"Permission denied: {file_path}")
                continue
            fp = os.path.abspath(file_path)
            file_md5_records.append((fp,md5_hash.hexdigest(), f"rm -rf {fp}"))
    df = pd.DataFrame(file_md5_records)
    df.columns = ["file", "md5", "cmd"]
    return df
@add_log
def run_task(dir = ""):
    out_name = dir.strip("/").split("/")[-1]
    run_task.logger.info(f"RUN {dir}")
    df = generate_md5_for_files(dir)
    df = df[df.duplicated(subset="md5", keep = False)]
    df = df.sort_values('md5')
    df.to_csv(f"./{out_name}_MD5.txt", sep = "\t", index = None)

@add_log
def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument(f'--path',help='', required=True) # 重复文件查找，每一行是一个目录
    args = parser.parse_args()    
    
    with open(f"{args.path}") as fd:
        tasks = [line.strip("\n") for line in fd.readlines()]
    
    #for t in tasks:
    #    run_task(t)
    
    main.logger.error(f"RUN {tasks}")
    results = []
    with ProcessPoolExecutor(16) as executor:
        for inx, task in enumerate(tasks):
            results.append(executor.submit(run_task, dir = task))
        for inx, res in enumerate(results): 
            try:
                res.result()
            except:
                main.logger.error(f"Error {tasks[inx]}")
    executor.shutdown(True)
    

if __name__ == "__main__":
    main()