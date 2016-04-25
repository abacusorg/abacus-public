# Ignore errors from the following command
# which is a hack to get install-sh
# (see http://www.freesoftwaremagazine.com/articles/configuring_a_project_with_autoconf)
automake --add-missing --copy &> /dev/null
autoreconf --install
