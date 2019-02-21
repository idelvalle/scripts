# NonacusTools

Go to folder containing NonacusTools:

```bash
cd ~/src/NonacusTools_1_Read
```

Look for folder with fastq files (named fastq) (include also batch_list.txt file)

Run the following command:

```bash
python consensus.py -i '/home/nacho/Desktop/Federica_Nonacus/fastq/batch_list.txt' -t 32
```

This will generate one folder for each sample.

Change permissions if necessary and move data:

```bash
sudo chown -R $USER A* # if first letter all folders is A
```

Call variants using platypus:

```bash
 bash PT.sh
#Be sure to have platypus.sh and PT.sh in same folder as files 
```

