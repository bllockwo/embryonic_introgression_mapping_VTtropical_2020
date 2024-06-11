### check file numbers

du -a | cut -d/ -f2 | sort | uniq -c | sort -nr