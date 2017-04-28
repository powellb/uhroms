for d in Compilers Data Lib Master ROMS User Waves; do
    echo "Replace ${d}"
    git rm -r ${d}
    cp -r ../trunk/${d} .
    git add ${d}
done
git commit -m "Update to latest ROMS trunk"    
