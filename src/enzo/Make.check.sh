#!/bin/bash
declare SCRIPT_NAME=$0
declare -i verbose=0
if (( $# == 1 ))
then
	read -a o_paths <<< $1
else
	declare -a o_paths=($@)
fi

# The outer product of mydirs and c_names is the list of
# source files to compile, minus the files to ignore, i_arr.
read -a mydirs <<< '. *'
read -a c_names <<< '*.c *.C *.F *.F90 *.cu *.pyx'
read -a c_paths <<< ''

# Specify more globs for files to be tracked by Mercurial in
# t_paths, in addition to c_names. The outer product of the 
# combined lists and dirs_globs is the list of files to track,
# minus the list of files to ignore, i_arr.
read -a t_names <<< '*.h *.cuh *.T *.py'
read -a t_paths <<< 'Makefile Make.*'

# The gnore-file list is made of two parts -- one that keeps the
# suffixes and another that has the suffixes removed.
declare -a i_names=( \
)
declare -a i_paths=( \
	'auto_show_*.C' \
	'auto_show_*.o' \
	'Make.config.machine' \
	'Make.config.override' \
)

if (( 0 )); then
	echo "Using test data......................"
	# data for testing without make
	read -a c_names <<< '*.cu *.pyx'
	read -a c_paths <<< 'arc*.C'
	read -a t_names <<< ''
	read -a t_paths <<< '*.T Makefile Make.* auto_*.C'
	read -a o_paths <<< " \
		arccosh.o \
		arcsinh.o \
		AssignGridToTaskMap.o \
		auto_show_compile_options.o \
	"
	i_names+=( \
		'no.such.file.*' \
		'*.pyx' \
	)
	#end test data
fi

pack(){
	declare -n arr2pack=$1
	arr2pack=(${arr2pack[*]})
}

expand_globs(){
	local -n a=$1

	# Perform elemntwise pathname expansion:
	pack a

	# Remove globs from the arra without a match:
	for i in ${!a[*]}; do 
		[ -e ${a[i]} ]
		(( $? == 1 )) && unset a[i]
	done

	# Pack indices
	pack a
}

del_suffixes(){
	local -n dest=$1
	local -n src=$2
	dest=(${src[@]%.*})
}

dirs_x_fnames(){
	local -n dest=$1
	local -n d_arr=$2
	local -n fn_arr=$3
	local -i i j
	local d

	if (( verbose )); then
		echo "$1 before:"
		echo "dirs: ${d_arr[@]}"
		echo "files: ${fn_arr[@]}"
	fi

	for i in ${!d_arr[@]}; do
		d="${d_arr[i]}/"

		if [ "$d" == "./" ]; then 
			d=""
		fi

		for j in ${!fn_arr[@]}; do
			dest+=("$d${fn_arr[j]}")
		done
	done

	if (( verbose )); then
		echo "$1 after: ${dest[@]}"
	fi
}

subtract_arrays(){
	local -n dest_arr=$1
	local -n l_arr=$2
	local -n r_arr=$3
	local -i i

	if [ "$1" != "$2" ]; then
		dest_arr=(${!l_arr[*]})
	fi

	# Concatenate arrays into strings:
	local r_str=" ${r_arr[*]} "

	if (( verbose )); then
		echo "$0: $2, left, ${#l_arr[*]}: ${l_arr[*]}"
		echo "$0: $3, right, ${#r_arr[*]}: $r_str"
	fi

	for i in ${!l_arr[*]}; do
		l_item=${l_arr[i]}
		if [[ $r_str == *" ${l_item} "* ]]; then
			unset dest_arr[i]
		fi
	done

	pack dest_arr

	if (( verbose )); then
		echo "$0: left/right, ${#dest_arr[*]}: ${dest_arr[*]}"
	fi

	return ${#dest_arr[*]}
}

print_arr_diff(){
	local -n n=$1
	local -n l_set=$2
	local -n r_set=$3
	local -n l_values=$4
	local -a dest_set=()
	local -i i j
	local print_cmd=$5
	local dest_str="-at"

	if (( verbose )); then
		echo "$0: $2=${l_set[*]}"
		echo "$0: $3=${r_set[*]}"
	fi

	subtract_arrays dest_set l_set r_set
	n=$?

	if (( n > 0 )); then
		echo ""
		if [ "$print_cmd" == "ls" ]; then
			echo "$n item(s) of total ${#l_set[@]} in $2 not in $3, newest first:"
			for i in ${!dest_set[*]}; do
				j=${dest_set[i]}
				dest_str+=" '${l_values[j]}'"
			done
			echo "$dest_str" | xargs ls -l
			echo ""
			echo "Same as a simple list, newest first:"
			dest_str=`echo "$dest_str" | xargs ls -1`
			dest_str=${dest_str/\n/}
			echo $dest_str
			
		else
			echo "$n item(s) of total ${#l_set[@]} in $2 not in $3:"
			for i in ${!dest_set[*]}; do
				j=${dest_set[i]}
				echo ${l_values[j]}
			done
		fi
	fi

	if (( verbose )); then
		echo "$0: $4(${#l_values[*]})=${l_values[*]}"
		echo "$0: !$4(${#l_values[*]})=${!l_values[*]}"
	fi

	return $n
}

echo "File patterns:"
echo "To compile, c_paths = (${c_paths[*]}) + (${mydirs[*]})*(${c_names[*]})"
echo "To ignore, i_paths  = (${i_paths[*]}) + (${mydirs[*]})*(${i_names[*]})"
echo "To track, t_paths   = (${t_paths[*]}) + (${mydirs[*]})*(${t_names[*]}) + c_paths"

dirs_x_fnames c_paths mydirs c_names
dirs_x_fnames i_paths mydirs i_names
dirs_x_fnames t_paths mydirs t_names

expand_globs c_paths
expand_globs i_paths
# o_paths shouldn't have any globs
expand_globs t_paths

if (( verbose )); then
	echo "After pathname expansion"
	echo "${#c_paths[*]} file(s) to compile"
	echo "${#o_paths[*]} object file(s)"
	echo "${#t_paths[*]} file(s) to track"
	echo "${#i_paths[*]} file(s) to ignore"
	#echo "c_paths=${c_paths[*]}"
	#echo "i_paths=${i_paths[*]}"
fi

subtract_arrays c_paths c_paths i_paths
subtract_arrays o_paths o_paths i_paths
subtract_arrays t_paths t_paths i_paths

t_paths+=(${c_paths[@]})

echo ""
echo "After pathname expansion and subtracting the files to ignore:"
echo "${#i_paths[*]} file(s) to ignore"
echo "${#c_paths[*]} file(s) to compile"
echo "${#o_paths[*]} object file(s)"
echo "${#t_paths[*]} file(s) to track"
echo "${#i_paths[*]} file(s) to have been ignored"
if (( verbose )); then
	echo "c_paths=${c_paths[*]}"
	echo "i_paths=${i_paths[*]}"
fi

## Remove suffixes:
declare -a c_paths2 o_paths2
del_suffixes c_paths2 c_paths
del_suffixes o_paths2 o_paths
 
declare -i co oc thg
print_arr_diff co c_paths2 o_paths2 c_paths 'ls'
print_arr_diff oc o_paths2 c_paths2 o_paths 'echo'
unset c_paths
unset o_paths

# Get paths of the files in $PWD tracked by hg.
# ( Note: hg outputs a list of paths relative to the repo root
#   directory.  Here it is converted to be relative to $PWD.  )
read -a hg_paths <<< `hg locate -I $PWD`
repo_root="`hg root`"
pwd_rel_to_root="${PWD#${repo_root}/}"
hg_paths=(${hg_paths[@]#${pwd_rel_to_root}/})
echo ""
echo "${#hg_paths[*]} files tracked by hg in '$rel_dir'."
if (( verbose )); then
	echo "The first 10:"
	for i in {1..10}; do
		echo ${hg_paths[i]}
	done
fi
print_arr_diff thg t_paths hg_paths t_paths 'ls'

echo ""
echo ""
echo "SUMMARY REPORT:"
echo ""
if (( $co )); then
	echo "$co exisiting source file(s) for compilation not matching names in"
	echo "  the objct files list.  Posibly working files need to be updated"
	echo "  to the correct revision. Try:"
	echo "      hg update tip"
	echo "      hg update -r <id>"
	echo "  where <id> is the desired revision id if not tip.  Also some"
	echo "  old files could had been excluded from the build."
else
	echo "All exisiting source files for compilation have a matching name in"
	echo "   the objct files list."
fi

echo ""
if (( $oc )); then
	echo "$oc object file name(s) not matching an existing source file."
	echo "  Consider editing $rel_dir/Make.config.objects and adding"
	echo "  the corresponding object file names to the proper list."
else
	echo "All object file names have a matching existing source file."
fi

echo ""
if (( thg )); then
	echo "$thg existing file(s) are not being tracked by hg but possibly"
	echo "  need to be.  See which files need to be added to the repo."
	echo "  Also make sure that the working revision is correct."
else
	echo "All existing files from the t_paths list are being tracked"
	echo "  by hg."
fi
echo "(End of summary)"
