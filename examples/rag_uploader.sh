#!/bin/bash

# RAG File Uploader for Google Cloud Vertex AI
# Taken from https://cloud.google.com/vertex-ai/generative-ai/docs/model-reference/rag-api#upload-a-rag-file-example-api

# Default values
PROJECT_ID="molviewAgent"
LOCATION="us-central1"
RAG_CORPUS_ID="4611686018427387904"
LOCAL_FILE_PATH=""
DISPLAY_NAME=""
DESCRIPTION=""

# Function to display usage
usage() {
    echo "Usage: $0 [OPTIONS]"
    echo ""
    echo "Upload a file to Google Cloud Vertex AI RAG corpus"
    echo ""
    echo "Required arguments:"
    echo "  -f, --file PATH          Local path to the file to be uploaded"
    echo "  -n, --name NAME          Display name of the RAG file"
    echo "  -d, --description DESC   Description of the RAG file"
    echo ""
    echo "Optional arguments:"
    echo "  -p, --project ID         Google Cloud Project ID (default: $PROJECT_ID)"
    echo "  -l, --location LOC       VertexAI location (default: $LOCATION)"
    echo "  -c, --corpus ID          RAG Corpus ID (default: $RAG_CORPUS_ID)"
    echo "  -h, --help               Show this help message"
    echo ""
    echo "Examples:"
    echo "  $0 -f document.pdf -n 'My Document' -d 'Important research paper'"
    echo "  $0 --file data.txt --name 'Dataset' --description 'Training data' --project myproject"
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -f|--file)
            LOCAL_FILE_PATH="$2"
            shift 2
            ;;
        -n|--name)
            DISPLAY_NAME="$2"
            shift 2
            ;;
        -d|--description)
            DESCRIPTION="$2"
            shift 2
            ;;
        -p|--project)
            PROJECT_ID="$2"
            shift 2
            ;;
        -l|--location)
            LOCATION="$2"
            shift 2
            ;;
        -c|--corpus)
            RAG_CORPUS_ID="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Error: Unknown option $1"
            usage
            exit 1
            ;;
    esac
done

# Validate required arguments
if [[ -z "$LOCAL_FILE_PATH" ]]; then
    echo "Error: File path is required (-f/--file)"
    usage
    exit 1
fi

if [[ -z "$DISPLAY_NAME" ]]; then
    echo "Error: Display name is required (-n/--name)"
    usage
    exit 1
fi

if [[ -z "$DESCRIPTION" ]]; then
    echo "Error: Description is required (-d/--description)"
    usage
    exit 1
fi

# Validate file exists
if [[ ! -f "$LOCAL_FILE_PATH" ]]; then
    echo "Error: File '$LOCAL_FILE_PATH' does not exist"
    exit 1
fi

# Check if gcloud is installed and authenticated
if ! command -v gcloud &> /dev/null; then
    echo "Error: gcloud CLI is not installed"
    exit 1
fi

# Check authentication
if ! gcloud auth print-access-token &> /dev/null; then
    echo "Error: Not authenticated with gcloud. Run 'gcloud auth application-default login' first"
    exit 1
fi

echo "Uploading file: $LOCAL_FILE_PATH"
echo "Display name: $DISPLAY_NAME"
echo "Description: $DESCRIPTION"
echo "Project: $PROJECT_ID"
echo "Location: $LOCATION"
echo "RAG Corpus ID: $RAG_CORPUS_ID"
echo ""

# Upload the file
curl -X POST \
    -H "X-Goog-Upload-Protocol: multipart" \
    -H "Authorization: Bearer $(gcloud auth print-access-token)" \
    -F metadata="{\"rag_file\": {\"display_name\":\"$DISPLAY_NAME\", \"description\":\"$DESCRIPTION\"}}" \
    -F file=@"$LOCAL_FILE_PATH" \
    "https://$LOCATION-aiplatform.googleapis.com/upload/v1beta1/projects/$PROJECT_ID/locations/$LOCATION/ragCorpora/$RAG_CORPUS_ID/ragFiles:upload"

echo ""
echo "Upload completed."
