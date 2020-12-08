sudo docker kill virus

sudo docker rm virus

sudo docker run --name virus -d -p 8383:5050 -v $PWD:/srv/shiny-server/ bioxfu/shiny-server

sudo docker exec -ti virus /bin/bash

# in docker container
# cd /srv/shiny-server/
# R
# shiny::runApp(host='0.0.0.0', port=5050, launch.browser=FALSE)
