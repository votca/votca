#
# votca-csg.spec
#
# Originally written by Jussi Lehtola <jussilehtola@fedoraproject.org>
# Fixed for multi-distro build by Klaus Kaempf <kkaempf@suse.de>
#
# Licensed under the Apache Software License (ASL 2.0)
#

Name:		votca-csg
Version:	1.0.1
Release:	1%{?dist}
Summary:	VOTCA coarse-graining engine
%if 0%{?suse_version}
Group:          Productivity/Scientific/Chemistry
%else
Group:		Applications/Engineering
%endif
License:	ASL 2.0
URL:		http://www.votca.org
Source0:	http://votca.googlecode.com/files/%{name}-%{version}.tar.bz2
Source1:	http://votca.googlecode.com/files/votca-tutorials.tar.bz2
Patch:          votca-csg-1.0.1.dif
BuildRoot:	%{_tmppath}/%{name}-%{version}-%{release}-root-%(%{__id_u} -n)

BuildRequires:	gromacs-devel
BuildRequires:  gcc-c++ libstdc++-devel
%if 0%{?suse_version}
BuildRequires:  pkg-config
BuildRequires:	libexpat-devel
%else
BuildRequires:  pkgconfig
BuildRequires:	expat-devel
%endif
BuildRequires:	fftw-devel
BuildRequires:	gsl-devel
BuildRequires:	boost-devel
BuildRequires:	libvotca-tools0-devel = %{version}

Requires:	%{name}-common = %{version}-%{release}
Requires:	%{name}-libs = %{version}-%{release}

%description
Versatile Object-oriented Toolkit for Coarse-graining Applications (VOTCA) is
a package intended to reduce the amount of routine work when doing systematic
coarse-graining of various systems. The core is written in C++.

This package contains the Coarse Graining Engine of VOTCA package.

%package libs
Summary:	Libraries for VOTCA coarse graining engine
%if 0%{?suse_version}
Group:		Productivity/Scientific/Chemistry
%else
Group:		System Environment/Libraries
%endif

%description libs
Versatile Object-oriented Toolkit for Coarse-graining Applications (VOTCA) is
a package intended to reduce the amount of routine work when doing systematic
coarse-graining of various systems. The core is written in C++.

This package contains libraries for the Coarse Graining Engine of VOTCA package.

%package devel
Summary:	Development headers and libraries for VOTCA Coarse Graining Engine
%if 0%{?suse_version}
Group:		Development/Libraries/C and C++
%else
Group:		Development/Libraries
%endif
Requires:	%{name}-libs = %{version}-%{release}
Requires:	votca-tools-devel = %{version}

%description devel
This package contains development headers and libraries for the Coarse Graining
Engine of VOTCA.

%package common
Summary:	Architecture independent data files for VOTCA CSG
%if 0%{?suse_version}
Group:          Productivity/Scientific/Chemistry
%else
Group:		Applications/Engineering
%endif
BuildArch:	noarch

%description common
Versatile Object-oriented Toolkit for Coarse-graining Applications (VOTCA) is
a package intended to reduce the amount of routine work when doing systematic
coarse-graining of various systems. The core is written in C++.

This package contains architecture independent data files for VOTCA CSG.

%package tutorials
Summary:	Tutorial documentation for VOTCA Coarse Graining Engine
%if 0%{?suse_version}
Group:          Productivity/Scientific/Chemistry
%else
Group:		Applications/Engineering
%endif
BuildArch:	noarch

%description tutorials
Versatile Object-oriented Toolkit for Coarse-graining Applications (VOTCA) is
a package intended to reduce the amount of routine work when doing systematic
coarse-graining of various systems. The core is written in C++.

This package contains tutorial documentation and sample data

%package bash
Summary:	Bash completion for votca
%if 0%{?suse_version}
Group:          Productivity/Other
%else
Group:		System Environment/Shells
%endif
Requires:	%{name} = %{version}-%{release}
Requires:	bash-completion
BuildArch:	noarch

%description bash
Versatile Object-oriented Toolkit for Coarse-graining Applications (VOTCA) is
a package intended to reduce the amount of routine work when doing systematic
coarse-graining of various systems. The core is written in C++. Iterative
methods are implemented using bash + perl.

This package contains bash completion support for votca-csg.

%prep
%setup -q
tar -xjf %{S:1}
rm -rf tutorials/.hg*
%patch -p1
autoreconf -f -i

%build
%if 0%{?suse_version}
%configure --disable-static --disable-la-files --disable-rc-files
%else
%configure --disable-static --disable-la-files --disable-rc-files
%endif
make %{?_smp_mflags}

%install
rm -rf %{buildroot}
make install DESTDIR=%{buildroot}
# Move bash completion file to correct location
mkdir -p %{buildroot}%{_sysconfdir}/bash_completion.d
cp scripts/csg-completion.bash %{buildroot}%{_sysconfdir}/bash_completion.d/votca

%if 0%{?suse_version}
%define pkgdocdir %{_docdir}/%{name}
%else
%define pkgdocdir %{_docdir}/%{name}-%{version}
%endif

%clean
rm -rf %{buildroot}

%post libs -p /sbin/ldconfig
%postun libs -p /sbin/ldconfig

%files
%defattr(-,root,root,-)
%doc CHANGELOG NOTICE README LICENSE
%{_bindir}/csg_*
%{_bindir}/multi_g_*

%files tutorials
%defattr(-,root,root,-)
%doc tutorials

%files common
%defattr(-,root,root,-)
%{_datadir}/votca

%files libs
%defattr(-,root,root,-)
%doc LICENSE
%{_libdir}/libvotca_csg.so.*

%files devel
%defattr(-,root,root,-)
%{_includedir}/votca/csg/
%{_libdir}/libvotca_csg.so
%{_libdir}/pkgconfig/libvotca_csg.pc

%files bash
%defattr(-,root,root,-)
%{_sysconfdir}/bash_completion.d/votca

%changelog
* Thu Nov 30 2010 Christoph Junghans <junghans@votca.org> - 1.0.1-1
- Minor cleanup.
* Thu Nov 25 2010 Jussi Lehtola <jussilehtola@fedoraproject.org> - 1.0-1
- First release.
