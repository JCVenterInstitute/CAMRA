#!/usr/bin/bash
#install if not alreadin installed
user_name="$1"
password="$2"
./login-expect "$user_name" "$password"