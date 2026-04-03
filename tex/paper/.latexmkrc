# Ensure the build directory exists
system("mkdir -p ./build");

# PDF generation mode (1 = pdflatex)
$pdf_mode = 1;

# Put all output (including PDF) in build directory first
$out_dir = 'build';

# Emulate aux directory behavior
$emulate_aux = 1;

# Clean up auxiliary files
$clean_ext = 'aux bbl blg fdb_latexmk fls log out synctex.gz toc';

# Copy PDF back to main directory after successful build
$success_cmd = 'cp build/%R.pdf ./';
