declare -A p1_arr
p1_arr["polblogs_d"]="25;0.9;0.9"
p1_arr["human_stelzl_d"]="10;0.52;0.54"
p1_arr["bitcoin_d"]="10;0.7;0.6"
p1_arr["mathoverflow_d"]="20;0.87;0.85"
p1_arr["epinions_d"]="20;0.8;0.9"
p1_arr["twitter_d"]="7;0.5;0.5"
p1_arr["web_google_d"]="15;0.75;0.8"
p1_arr["baidu_d"]="20;0.7;0.8"
p1_arr["flickr_d"]="9999;0.999;0.9999"
p1_arr["dbpedia_d"]="9999;0.999;0.9999"
p1_arr["usa_road_d"]="7;0.5;0.5"
p1_arr["clue_web_d"]="35;0.9;0.9"

compers="32"
tau_split="200"
tau_time="1"

datasets="polblogs_d human_stelzl_d bitcoin_d mathoverflow_d epinions_d twitter_d web_google_d baidu_d flickr_d dbpedia_d usa_road_d clue_web_d"
for f in $datasets
do

	mkdir -p ${f}

	p1_value=${p1_arr[$f]} # get p1 values by dataset name
	arr_split=(${p1_value//;/ }) # split values by semi colon ;
	t_size=${arr_split[0]}
	gamma1=${arr_split[1]}
	gamma2=${arr_split[2]}

	base=${f}/"${compers}"
	log="${base}_log"
	(date;time ./run ${f} ${compers} ${gamma1} ${gamma2} ${t_size} ${tau_time} ${tau_split};date) 2>&1 | tee ${log};
	./run_scala.sh ${d} ${c}
done
