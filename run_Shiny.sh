docker kill virus

docker rm virus

docker run --name virus -d -p 8383:3838 -v $PWD:/srv/shiny-server/ bioxfu/shiny-server
