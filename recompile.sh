PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC_DIR="${PROJECT_DIR}/src"

cd "${SRC_DIR}" || exit 1
make clean
make
cd ..