#add colors to echo's
YELLOW='\033[1;33m'
GREEN='\033[1;32m'
RED='\033[1;31m'
RESET='\033[0m'

# types of sintax:
# ./run.sh [nu_1] [nu_2] [plot_type]
# ./run.sh


# if no arguments are given, the default values are used (nu_1 = 0.5, nu_2 = 0.5, plot_type = 1)
if [ $# -eq 0 ]; then
  echo -e "${YELLOW}No arguments given! Using default values.${RESET}"
  nu_1=0.5
  nu_2=0.5
  plot_type=1
else
  nu_1=$1
  nu_2=$2
  plot_type=$3
fi

# if nu_1 <= 0.0 or nu_2 <= 0.0 or plot_type is not 1, 2 or 3, the program exits
if [ $nu_1 \< 0.0 ] || [ $nu_2 \< 0.0 ] || ([ $plot_type != 1 ] && [ $plot_type != 2 ] && [ $plot_type != 3 ]); then
  echo -e "${RED}Invalid argument! Skipped animating.${RESET}"
  anim=0
else
  anim=1
fi

echo -e "${YELLOW}Compiling...${RESET}"
# make
make
if [ $? -ne 0 ]; then
  echo -e "${RED}Compilation failed!${RESET}"
  exit 1
fi
echo -e "${GREEN}Compiling done!${RESET}"
echo -e "${YELLOW}Running...${RESET}"
./bin/main $nu_1 $nu_2
if [ $? -ne 0 ]; then
  echo -e "${RED}Running failed!${RESET}"
  exit 1
fi
echo -e "${GREEN}Running done!${RESET}"

# do not plot if there is an extra argument to ./run.sh
if [ $anim -eq 1 ]; then
  echo -e "${YELLOW}Plotting...${RESET}"
  python src/plot.py $nu_1 $nu_2 $plot_type
  if [ $? -ne 0 ]; then
    echo -e "${RED}Plotting failed!${RESET}"
    exit 1
  fi
  echo -e "${GREEN}Plotting done!${RESET}"
fi

echo -e "${YELLOW}Plotting energy...${RESET}"
python src/plot_energy.py $nu_1 $nu_2
if [ $? -ne 0 ]; then
  echo -e "${RED}Plotting energy failed!${RESET}"
  exit 1
fi
