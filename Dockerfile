# The `autorift_itslive` bundle on a bare Debian: the bundle carries every library it loads but
# the C library. The build context holds each platform's unpacked bundle as `<arch>/autorift_itslive`
# (`amd64`, `arm64`, as Docker names them), so one build makes either image, or both:
#
#   docker buildx build --platform linux/amd64,linux/arm64 -t autorift_itslive .
#
# The release workflow (`.github/workflows/Release.yml`) stages the bundles and builds it.
FROM debian:bookworm-slim

ARG TARGETARCH
COPY ${TARGETARCH}/autorift_itslive /opt/autorift_itslive
# A home of its own, writable by the user the image runs as: a product is written to the working directory
# when no output path is given.
RUN useradd --uid 1000 --create-home autorift

ENV PATH=/opt/autorift_itslive/bin:$PATH
USER 1000
WORKDIR /home/autorift

ENTRYPOINT ["autorift_itslive"]
CMD ["--help"]
