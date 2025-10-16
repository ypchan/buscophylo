ls scripts | grep '.py' | while read s; do
    rm ~/bin/${s}
    chomd 755 scripts/${s}
    ln -s $(pwd)/scripts/${s} ~/bin/${s}
    echo "ln -s $(pwd)/scripts/${s} ~/bin/${s}"
done