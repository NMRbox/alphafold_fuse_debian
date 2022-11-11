# Debian package

dpkg-buildpackage -b -uc -us builds the package.

This includes ensuring libfuse-dev is installed.

The mount directory is created on package installation by debian/dirs.

debian/postinst and debian/prerm install and remove the alphafoldfuse service, respectively.
