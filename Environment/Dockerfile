# Start from the official Perl base image
FROM perl:latest

# Install cpanm (if not included in the image already)
RUN curl -L https://cpanmin.us | perl - --self-upgrade

# Download and install Diamond BLASTX
RUN wget http://github.com/bbuchfink/diamond/releases/download/v2.1.10/diamond-linux64.tar.gz \
    && tar -xzf diamond-linux64.tar.gz \
    && mv diamond /usr/local/bin/ \
    && rm diamond-linux64.tar.gz 

# Set up the working directory (optional)
WORKDIR /app

# Copy the scripts from the local machine to the Docker image
COPY ./create-labeler.pl /app
COPY ./kraken-parser.pl /app

# Make sure the scripts are executable
RUN chmod +x /app/create-labeler.pl /app/kraken-parser.pl

# Install Perl modules using cpanm (optional)
RUN cpanm --force Bio::MUST::Apps::Physeter

# Add /app to PATH
ENV PATH=$PATH:/app

# Set the default command (optional)
CMD ["perl", "-v"]