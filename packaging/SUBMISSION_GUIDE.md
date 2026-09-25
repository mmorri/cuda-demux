# Package Submission Guide for cuda-demux

This guide explains how to submit cuda-demux packages to official distribution repositories.

## 📦 Arch Linux (AUR)

### Prerequisites
- AUR account at https://aur.archlinux.org/register
- SSH key added to your AUR account
- `git` and `makepkg` installed

### Submission Steps

1. **Test the package locally:**
```bash
cd packaging/arch/
makepkg -si
# Test that cuda-demux works
cuda-demux --help
```

2. **Create AUR repository:**
```bash
# Clone empty AUR repository
git clone ssh://aur@aur.archlinux.org/cuda-demux.git aur-cuda-demux
cd aur-cuda-demux
```

3. **Copy and update PKGBUILD:**
```bash
cp ../packaging/arch/PKGBUILD .
# Generate .SRCINFO
makepkg --printsrcinfo > .SRCINFO
```

4. **Submit to AUR:**
```bash
git add PKGBUILD .SRCINFO
git commit -m "Initial submission of cuda-demux 1.0.0"
git push
```

5. **Maintain the package:**
- Monitor comments at https://aur.archlinux.org/packages/cuda-demux
- Update for new releases
- Respond to user issues

### AUR Guidelines
- Follow: https://wiki.archlinux.org/title/AUR_submission_guidelines
- Use namcap to check package: `namcap PKGBUILD`
- Test in clean chroot: `extra-x86_64-build`

---

## 🐧 Debian/Ubuntu

### Official Debian Repository

1. **Become a Debian Maintainer:**
   - Guide: https://mentors.debian.net/intro-maintainers
   - Need a Debian Developer sponsor

2. **Prepare the package:**
```bash
# Install packaging tools
sudo apt install devscripts dh-make lintian

# Create source package
cd cuda-demux-1.0.0/
dpkg-buildpackage -S -sa

# Check package quality
lintian ../cuda-demux_1.0.0-1_amd64.changes
```

3. **Submit to mentors.debian.net:**
   - Create account at https://mentors.debian.net/
   - Upload package using dput:
```bash
dput mentors ../cuda-demux_1.0.0-1_amd64.changes
```

4. **Request sponsorship:**
   - File RFS (Request For Sponsorship) bug
   - Email debian-mentors@lists.debian.org

### Ubuntu PPA (Personal Package Archive)

**Easier alternative for Ubuntu users:**

1. **Create Launchpad account:**
   - Register at https://launchpad.net/
   - Add GPG key and SSH key

2. **Create PPA:**
```bash
# On Launchpad web interface:
# Click "Create a new PPA"
# Name: cuda-demux
```

3. **Build and upload:**
```bash
# Update changelog for Ubuntu
dch -D focal -v 1.0.0-1ubuntu1 "Initial Ubuntu release"

# Build source package
debuild -S -sa

# Sign and upload to PPA
dput ppa:yourusername/cuda-demux ../cuda-demux_1.0.0-1ubuntu1_source.changes
```

4. **Users can install from PPA:**
```bash
sudo add-apt-repository ppa:yourusername/cuda-demux
sudo apt update
sudo apt install cuda-demux
```

---

## 🎩 Fedora/RHEL

### Fedora Official Repository

1. **Become a Fedora Package Maintainer:**
   - Create FAS account: https://accounts.fedoraproject.org/
   - Join packager group: https://docs.fedoraproject.org/en-US/package-maintainers/Joining_the_Package_Maintainers/

2. **Submit package for review:**
```bash
# Create SRPM
rpmbuild -bs packaging/rpm/cuda-demux.spec

# Upload SRPM to fedorapeople.org
scp ~/rpmbuild/SRPMS/cuda-demux-1.0.0-1.src.rpm \
    yourusername@fedorapeople.org:~/public_html/
```

3. **File review request:**
   - Go to https://bugzilla.redhat.com/
   - Product: Fedora
   - Component: Package Review
   - Title: "Review Request: cuda-demux - GPU-accelerated Illumina demultiplexer"
   - Include SPEC and SRPM URLs

4. **Pass review process:**
   - Address reviewer feedback
   - Run fedora-review tool
   - Get approved by sponsor

### Fedora COPR (Easier Alternative)

**Community repository - easier to set up:**

1. **Create COPR account:**
   - Login at https://copr.fedorainfracloud.org/
   - Create new project

2. **Build from Git:**
```bash
# In COPR web interface:
# Packages → New Package
# Package name: cuda-demux
# Clone URL: https://github.com/mmorri/cuda-demux.git
# Spec file: packaging/rpm/cuda-demux.spec
```

3. **Users can install:**
```bash
sudo dnf copr enable yourusername/cuda-demux
sudo dnf install cuda-demux
```

---

## 🔧 Pre-submission Checklist

### All Distributions:
- [ ] Package builds successfully
- [ ] Binary runs and produces correct output
- [ ] Dependencies are correctly specified
- [ ] License file is included
- [ ] Documentation is included
- [ ] Version numbering follows distribution conventions

### Arch Linux:
- [ ] PKGBUILD passes `namcap` checks
- [ ] .SRCINFO is generated and included
- [ ] Package name doesn't conflict with existing packages

### Debian/Ubuntu:
- [ ] Package passes `lintian` checks
- [ ] debian/copyright file is complete
- [ ] debian/watch file for upstream monitoring
- [ ] Follows Debian Policy Manual

### Fedora/RHEL:
- [ ] Spec file passes `rpmlint` checks
- [ ] Follows Fedora Packaging Guidelines
- [ ] Source URL is accessible
- [ ] BuildRequires are complete

---

## 📝 Package Metadata Updates

When submitting, update these fields in the package files:

1. **Maintainer Information:**
   - Replace "Your Name <your.email@example.com>" with actual details

2. **Version Tags:**
   - Create Git tags for releases: `git tag v1.0.0`
   - Update version in all package files

3. **Checksums (Arch only):**
   - Generate: `makepkg -g`
   - Update sha256sums in PKGBUILD

---

## 🚀 Quick Start Commands

### Build all packages locally:
```bash
cd packaging/
./build-packages.sh all
```

### Test installations:
```bash
# Arch
sudo pacman -U arch/*.pkg.tar.zst

# Debian/Ubuntu
sudo dpkg -i *.deb

# Fedora/RHEL
sudo rpm -ivh *.rpm
```

---

## 📚 Resources

### Arch Linux
- AUR Guidelines: https://wiki.archlinux.org/title/AUR_submission_guidelines
- PKGBUILD Reference: https://wiki.archlinux.org/title/PKGBUILD

### Debian
- New Maintainer's Guide: https://www.debian.org/doc/manuals/maint-guide/
- Debian Policy: https://www.debian.org/doc/debian-policy/
- Mentors FAQ: https://mentors.debian.net/intro-maintainers

### Ubuntu
- Packaging Guide: https://packaging.ubuntu.com/
- PPA Help: https://help.launchpad.net/Packaging/PPA

### Fedora
- Package Maintainer Docs: https://docs.fedoraproject.org/en-US/package-maintainers/
- Packaging Guidelines: https://docs.fedoraproject.org/en-US/packaging-guidelines/
- COPR Documentation: https://docs.pagure.org/copr.copr/

---

## ⚠️ Important Notes

1. **CUDA Dependencies:** Some distributions may not have CUDA in official repos. Consider:
   - Documenting NVIDIA's CUDA repository setup
   - Using conditional dependencies
   - Providing static builds as alternative

2. **Testing:** Always test in clean environments:
   - Docker containers
   - Virtual machines
   - Mock/sbuild/pbuilder

3. **Licensing:** Ensure all dependencies are compatible with MIT license

4. **Security:** Enable reproducible builds where possible

5. **Updates:** Set up CI/CD for automatic package building on new releases