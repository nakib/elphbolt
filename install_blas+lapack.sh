#!/bin/bash

cd thirdparty

cd forblas
fpm install --profile release


cd ../forlapack
fpm install --profile release

#This will create the libraries libforblas.a and libforlapack.a in ~/.local/lib/
