#!/usr/bin/expect -f

# Get the arguments passed from the bash script
set user_name [lindex $argv 0]
set password [lindex $argv 1]

set timeout -1

# Start the login process
spawn p3-login $user_name
expect "password: "
send -- "$password\r"

# Allow user interaction
interact