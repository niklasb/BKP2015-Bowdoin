#!/bin/bash
grep -v -- ----- | tr -d '\n' | base64 -d | openssl asn1parse -inform DER -i
