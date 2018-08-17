sudo docker kill virus

sudo docker rm virus

sudo docker run --name virus -d -p 8383:3838 -v $PWD:/srv/shiny-server/ bioxfu/shiny-server
