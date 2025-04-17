#!/bin/bash

# 编译a.cpp到a.o
# g++ a.cpp -o a.o
g++ a.cpp -o a.o -DDEBUG
if [ $? -ne 0 ]; then
    echo "Compilation Failed"
    exit 1
fi

# 运行程序并比较输出
./a.o < a.in > output.tmp

echo "finished."

# 比较输出文件与预期结果
if cmp -s output.tmp a.out; then
    echo "Accepted"
else
    echo "Wrong Answer"
#     echo ""
#     echo "=== Expected Output ==="
#     cat a.out
#     echo ""
#     echo "=== Your Output ==="
#     cat output.tmp
#     echo ""
fi

# 清理临时文件
# rm -f output.tmp