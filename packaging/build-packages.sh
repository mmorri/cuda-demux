#!/bin/bash
# Build script for creating distribution packages

set -e

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PROJECT_ROOT="$( cd "$SCRIPT_DIR/.." && pwd )"
VERSION="1.0.0"

echo "Building cuda-demux packages version $VERSION"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Function to print colored output
print_status() {
    echo -e "${GREEN}[*]${NC} $1"
}

print_error() {
    echo -e "${RED}[!]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[!]${NC} $1"
}

# Build AUR package
build_aur() {
    print_status "Building AUR package..."
    cd "$SCRIPT_DIR/arch"

    # Update version in PKGBUILD
    sed -i "s/pkgver=.*/pkgver=$VERSION/" PKGBUILD

    # Create source tarball
    cd "$PROJECT_ROOT"
    git archive --format=tar.gz --prefix=cuda-demux-$VERSION/ HEAD > "$SCRIPT_DIR/arch/cuda-demux-$VERSION.tar.gz"

    # Build package
    cd "$SCRIPT_DIR/arch"
    makepkg -sf --noconfirm || print_warning "AUR build failed - you may need to install dependencies"

    print_status "AUR package built: $(ls *.pkg.tar.zst 2>/dev/null || echo 'Not found')"
}

# Build Debian package
build_deb() {
    print_status "Building Debian package..."

    # Create temporary build directory
    BUILD_DIR="/tmp/cuda-demux-$VERSION"
    rm -rf "$BUILD_DIR"
    cp -r "$PROJECT_ROOT" "$BUILD_DIR"

    # Copy debian directory
    cp -r "$SCRIPT_DIR/debian" "$BUILD_DIR/"

    # Update changelog version
    cd "$BUILD_DIR"
    sed -i "s/cuda-demux (.*/cuda-demux ($VERSION-1) unstable; urgency=medium/" debian/changelog

    # Build package
    dpkg-buildpackage -us -uc -b || print_warning "DEB build failed - you may need to install build dependencies"

    # Copy packages back
    cp /tmp/cuda-demux*.deb "$SCRIPT_DIR/" 2>/dev/null || true

    print_status "DEB package built: $(ls $SCRIPT_DIR/*.deb 2>/dev/null || echo 'Not found')"
}

# Build RPM package
build_rpm() {
    print_status "Building RPM package..."

    # Setup RPM build tree
    mkdir -p ~/rpmbuild/{BUILD,RPMS,SOURCES,SPECS,SRPMS}

    # Create source tarball
    cd "$PROJECT_ROOT"
    git archive --format=tar.gz --prefix=cuda-demux-$VERSION/ HEAD > ~/rpmbuild/SOURCES/cuda-demux-$VERSION.tar.gz

    # Copy spec file
    cp "$SCRIPT_DIR/rpm/cuda-demux.spec" ~/rpmbuild/SPECS/

    # Update version in spec
    sed -i "s/Version:.*/Version:        $VERSION/" ~/rpmbuild/SPECS/cuda-demux.spec

    # Build RPM
    cd ~/rpmbuild
    rpmbuild -ba SPECS/cuda-demux.spec || print_warning "RPM build failed - you may need to install build dependencies"

    # Copy packages back
    cp ~/rpmbuild/RPMS/x86_64/cuda-demux*.rpm "$SCRIPT_DIR/" 2>/dev/null || true
    cp ~/rpmbuild/SRPMS/cuda-demux*.src.rpm "$SCRIPT_DIR/" 2>/dev/null || true

    print_status "RPM package built: $(ls $SCRIPT_DIR/*.rpm 2>/dev/null || echo 'Not found')"
}

# Main script
case "${1:-all}" in
    arch|aur)
        build_aur
        ;;
    deb|debian)
        build_deb
        ;;
    rpm|redhat|fedora)
        build_rpm
        ;;
    all)
        build_aur
        build_deb
        build_rpm
        ;;
    *)
        echo "Usage: $0 [arch|deb|rpm|all]"
        exit 1
        ;;
esac

print_status "Package building complete!"
print_status "Packages can be found in: $SCRIPT_DIR"