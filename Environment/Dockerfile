# Start from the official Perl base image
FROM perl:latest

# Install cpanm (if not included in the image already)
RUN curl -L https://cpanmin.us | perl - --self-upgrade

# Set up the working directory (optional)
WORKDIR /app

# Install Perl modules using cpanm (optional)
RUN cpanm --force Bio::MUST::Apps::Physeter

# Set the default command (optional)
CMD ["perl", "-v"]