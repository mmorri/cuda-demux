Name:           cuda-demux
Version:        1.0.0
Release:        1%{?dist}
Summary:        GPU-accelerated Illumina sequencing demultiplexer

License:        MIT
URL:            https://github.com/mmorri/cuda-demux
Source0:        https://github.com/mmorri/cuda-demux/archive/v%{version}/%{name}-%{version}.tar.gz

BuildRequires:  cmake >= 3.16
BuildRequires:  gcc-c++ >= 7
BuildRequires:  cuda-toolkit-11-0 >= 11.0
BuildRequires:  zlib-devel
BuildRequires:  libomp-devel
BuildRequires:  git
BuildRequires:  pkgconfig

Requires:       cuda-runtime-11-0 >= 11.0
Requires:       cuda-cudart-11-0 >= 11.0
Requires:       zlib >= 1.2.11
Requires:       libomp >= 5.0
Requires:       glibc >= 2.34
Requires:       libstdc++ >= 11
Requires:       libgcc

%description
CUDA-Demux is a high-performance, GPU-accelerated tool for demultiplexing
Illumina sequencing data from BCL/CBCL format to FASTQ files. It leverages
CUDA for fast barcode matching and parallel processing.

Features:
- GPU Acceleration: Fast barcode matching using CUDA kernels
- CBCL Support: Native support for compressed BCL files
- Multi-threaded: Parallel BCL/CBCL file parsing and decompression
- Dual Indexing: Support for single and dual-indexed libraries
- High Performance: Process billions of reads efficiently

%prep
%autosetup -n %{name}-%{version}

%build
mkdir -p build
cd build
cmake .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX=%{_prefix} \
    -DENABLE_FETCH_TINYXML2=ON

%make_build

%install
cd build
install -Dm755 cuda-demux %{buildroot}%{_bindir}/cuda-demux
cd ..
install -Dm644 README.md %{buildroot}%{_docdir}/%{name}/README.md
install -Dm644 LICENSE %{buildroot}%{_licensedir}/%{name}/LICENSE

%files
%license LICENSE
%doc README.md
%{_bindir}/cuda-demux

%changelog
* Thu Sep 25 2025 Your Name <your.email@example.com> - 1.0.0-1
- Initial package release
- GPU-accelerated BCL/CBCL to FASTQ demultiplexing
- Support for NovaSeq, NextSeq, MiSeq platforms
- Adaptive batch sizing for GPU memory management