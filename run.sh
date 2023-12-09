#add colors to echo's
YELLOW='\033[1;33m'
GREEN='\033[1;32m'
RED='\033[1;31m'
RESET='\033[0m'

echo -e "${YELLOW}Compiling...${RESET}"
# make
make
if [ $? -ne 0 ]; then
  echo -e "${RED}Compilation failed!${RESET}"
  exit 1
fi
echo -e "${GREEN}Compiling done!${RESET}"
echo -e "${YELLOW}Running...${RESET}"
./bin/main
if [ $? -ne 0 ]; then
  echo -e "${RED}Running failed!${RESET}"
  exit 1
fi
echo -e "${GREEN}Running done!${RESET}"

# do not plot if there is an extra argument to ./run.sh
if [ $# -eq 0 ]; then
  echo -e "${YELLOW}Plotting...${RESET}"
  python src/plot.py
  if [ $? -ne 0 ]; then
    echo -e "${RED}Plotting failed!${RESET}"
    exit 1
  fi
  echo -e "${GREEN}Plotting done!${RESET}"
else
  echo -e "${YELLOW}Skipped plotting.${RESET}"
fi