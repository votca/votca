#
# libvotca-tools0.spec
#
# Originally written by Jussi Lehtola <jussilehtola@fedoraproject.org>
# Fixed for multi-distro build by Klaus Kaempf <kkaempf@suse.de>
#
# Licensed under the Apache Software License (ASL 2.0)
#

Name:		libvotca-tools0
%define srcname votca-tools
Version:	1.0.1
Release:	1%{?dist}
Summary:	VOTCA tools library
%if 0%{?suse_version}
Group:          Productivity/Scientific/Chemistry
%else
Group:		Applications/Engineering
%endif
License:	ASL 2.0
URL:		http://www.votca.org
Source0:	http://votca.googlecode.com/files/%{srcname}-%{version}.tar.bz2

BuildRoot:	%{_tmppath}/%{name}-%{version}-%{release}-root-%(%{__id_u} -n)

BuildRequires:  gcc-c++ libstdc++-devel
%if 0%{?suse_version}
BuildRequires:  pkg-config
BuildRequires:	libexpat-devel
BuildRequires:	fftw3-devel
%else
BuildRequires:  pkgconfig
BuildRequires:	fftw-devel
BuildRequires:	expat-devel
%endif
BuildRequires:	gsl-devel
BuildRequires:	boost-devel

%description
Versatile Object-oriented Toolkit for Coarse-graining Applications (VOTCA) is
a package intended to reduce the amount of routine work when doing systematic
coarse-graining of various systems. The core is written in C++.

This package contains the basic tools library of VOTCA package.

%package devel
Summary:	Development headers and libraries for votca-tools
%if 0%{?suse_version}
Group:          Development/Libraries/C and C++
%else
Group:		Development/Libraries
%endif
Requires:	%{name} = %{version}-%{release}

%description devel
Versatile Object-oriented Toolkit for Coarse-graining Applications (VOTCA) is
a package intended to reduce the amount of routine work when doing systematic
coarse-graining of various systems. The core is written in C++.

This package contains development headers and libraries for votca-tools.

%prep
%setup -q -n %{srcname}-%{version}
# Get rid of bundled versions of boost and expat
rm -rf src/libboost
rm -rf src/libexpat
autoreconf -f -i

%build
%configure --disable-static --disable-la-files --disable-rc-files
make %{?_smp_mflags}

%install
rm -rf %{buildroot}
make install DESTDIR=%{buildroot}

%clean
rm -rf %{buildroot}

%post -p /sbin/ldconfig
%postun -p /sbin/ldconfig

%files
%defattr(-,root,root,-)
%doc CHANGELOG LICENSE NOTICE
%{_libdir}/libvotca_tools.so.*

%files devel
%defattr(-,root,root,-)
%{_includedir}/votca/
%{_libdir}/libvotca_tools.so
%{_libdir}/pkgconfig/libvotca_tools.pc

%changelog
* Thu Nov 30 2010 Christoph Junghans <junghans@votca.org> - 1.0.1-1
- Minor cleanup.
* Thu Nov 25 2010 Jussi Lehtola <jussilehtola@fedoraproject.org> - 1.0-1
- First release.
