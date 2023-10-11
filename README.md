# Efficient Approximation Framework for Attribute Recommendation
## Input Table Format
The table file needs to follow this rule:

1. The first line contains column names, which are strings seperated by commas.
2. Other lines are records, where each line consists of integers seperated by commas.

Example Table file:
```
ID,Age,IncomePerHour
1,35,10
2,25,23
3,59,13
4,61,42
5,27,33
6,30,45
```
A table with 3 attributes and 6 records.
## Example Usage
```
workSpace="yourWorkSpace"

# build the code
${workSpace}/ make clean
${workSpace}/ make

# serialize the original datafile
${workSpace}/ ./test save_serialize_file --prefix /data/ --dataset hus0920

# run a query
${workSpace}/ ./test topk_query --metric ks --prefix /data/ --dataset hus0920 --abs_error 0.05 --topk 4 --query_num 10 --initial_size 1024
```
## Parameters
--prefix &lt;prefix&gt;

--dataset &lt;dataset&gt;

--support_file &lt;support_file&gt;

--abs_error &lt;error_bound&gt;

--topk &lt;topk&gt;

--initial_size &lt;intial_size&gt;

--query_num &lt;query_num&gt;

--metric &lt;metric_function&gt;

For more information about parameters, you can run 
```
${workSpace}/ ./test test --help
```
