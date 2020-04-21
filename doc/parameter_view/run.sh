#!/bin/bash

echo "opening http://localhost:8011/parameters.xml"
xdg-open "http://localhost:8011/parameters.xml" &
python3 -m http.server 8011
