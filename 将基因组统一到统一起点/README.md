通过基因组文件获取每个基因的长度然后得到一个字典，key是基因名，value是元组，每个基因的累加（start，start+length），每个基因的start是上一个基因end-1，累加通过numpy的cumsum实现

然后将基因匹配进去即可，将start和end分别和字典里的start和end相加即可